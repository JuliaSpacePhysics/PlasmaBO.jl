# Electromagnetic fluid dispersion relation solver
# Ported from MATLAB bo_fluid_matrix.m by Hua-sheng XIE
#
# References:
# - Xie, H. S. (2014). PDRF: A general dispersion relation solver for
#   magnetized multi-fluid plasma. Comput. Phys. Comm. 185, 670-675.

export FluidSpecies, FluidSolverParams, solve_fluid_dispersion

"""
    FluidSpecies{T}

Fluid species parameters for the multi-fluid dispersion relation solver.

# Fields
- `q`: Charge in units of electron charge
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity (in units of c)
- `gamma_z`: Parallel polytrope exponent (default: 1.0 for isothermal)
- `gamma_p`: Perpendicular polytrope exponent (default: 1.0 for isothermal)
"""
struct FluidSpecies{T}
    q::T           # Charge (units of e)
    m::T           # Mass (units of mp)
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    gamma_z::T     # Parallel polytrope exponent
    gamma_p::T     # Perpendicular polytrope exponent
end

charge(s) = s.q * q
mass(s) = s.m * mp

plasma_frequency(s::FluidSpecies) = plasma_frequency(charge(s), s.n, mass(s))
gyrofrequency(s::FluidSpecies) = gyrofrequency(charge(s), s.B, mass(s))

"""
    FluidSpecies(q, m, n, Tz, Tp; vdz=0.0, gamma_z=1.0, gamma_p=1.0)

Create a FluidSpecies with specified parameters.

# Arguments
- `q`: Charge in units of e (-1 for electron, 1 for proton)
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity in units of c (default: 0)
- `gamma_z`: Parallel polytrope exponent (default: 1.0, isothermal)
- `gamma_p`: Perpendicular polytrope exponent (default: 1.0, isothermal)

# Fluid Closure Models
- `gamma_z = gamma_p = 1.0`: Isothermal
- `gamma_z = 3.0, gamma_p = 2.0`: CGL (Chew-Goldberger-Low) model

See also: [`solve_fluid_dispersion`](@ref)
"""
function FluidSpecies(q, m, n, Tz = 0.0, Tp = 0.0; vdz = 0.0, gamma_z = 1.0, gamma_p = 1.0)
    T = promote_type(
        typeof(q), typeof(m), typeof(n), typeof(Tz), typeof(Tp),
        typeof(vdz), typeof(gamma_z), typeof(gamma_p)
    )
    return FluidSpecies{T}(T(q), T(m), T(n), T(Tz), T(Tp), T(vdz), T(gamma_z), T(gamma_p))
end

"""
    FluidSolverParams{T}

Parameters for the fluid dispersion relation solver.
"""
struct FluidSolverParams{T} <: AbstractSolverParams
    qs::Vector{T}       # Charges (C)
    ms::Vector{T}       # Masses (kg)
    ns0::Vector{T}      # Number densities (m^-3)
    vdsz::Vector{T}     # Parallel drift velocities (m/s)
    Psz::Vector{T}      # Parallel pressures (Pa)
    Psp::Vector{T}      # Perpendicular pressures (Pa)
    gamma_z::Vector{T}  # Parallel polytrope exponents
    gamma_p::Vector{T}  # Perpendicular polytrope exponents
    B0::T               # Magnetic field (T)
end

"""
    create_fluid_params(species::Vector{<:FluidSpecies}, B0)

Create FluidSolverParams from a list of FluidSpecies.

# Arguments
- `species`: Vector of FluidSpecies objects
- `B0`: Magnetic field strength (Tesla)
"""
function create_fluid_params(species, B0)
    T = Float64

    qs = charge.(species)
    ms = mass.(species)
    ns0 = [s.n for s in species]

    # Temperatures in Kelvin
    Tzs = [s.Tz * q / kb for s in species]
    Tps = [s.Tp * q / kb for s in species]

    # Drift velocities
    vdsz = [s.vdz * c0 for s in species]

    # Pressures
    Psz = @. kb * ns0 * Tzs
    Psp = @. kb * ns0 * Tps

    # Polytrope exponents
    gamma_z = [s.gamma_z for s in species]
    gamma_p = [s.gamma_p for s in species]

    return FluidSolverParams{T}(
        T.(qs), T.(ms), T.(ns0), T.(vdsz),
        T.(Psz), T.(Psp), T.(gamma_z), T.(gamma_p), T(B0)
    )
end

"""
    solve_fluid_dispersion(params::FluidSolverParams, kx, kz)

Solve the multi-fluid electromagnetic dispersion relation using matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

The fluid model includes:
- Continuity equation: ∂n/∂t + ∇·(nv) = 0
- Momentum equation: m(∂v/∂t + v·∇v) = q(E + v×B) - ∇P/n
- Maxwell's equations for E and B

# Arguments
- `params`: FluidSolverParams with plasma parameters
- `kx`: Perpendicular wave vector component (m⁻¹)
- `kz`: Parallel wave vector component (m⁻¹)

# Returns
- Vector of complex eigenfrequencies ω (rad/s)

See also: [`FluidSpecies`](@ref), [`create_fluid_params`](@ref)
"""
function solve_fluid_dispersion(params::FluidSolverParams, kx, kz)
    (; qs, ms, ns0, vdsz, Psz, Psp, gamma_z, gamma_p, B0) = params
    S = length(qs)
    # Matrix dimensions: 4 variables per species (n, vx, vy, vz) + 6 fields (Ex, Ey, Ez, Bx, By, Bz)
    SJ = 4 * S
    NN = SJ + 6

    # Build sparse matrix
    M = zeros(ComplexF64, NN, NN)

    for s in 1:S
        ind = (s - 1) * 4
        qⱼ = qs[s]
        mⱼ = ms[s]
        nⱼ = ns0[s]
        ρⱼ = mⱼ * nⱼ # Mass densities (kg/m^3)
        wcⱼ = qⱼ * B0 / mⱼ # Cyclotron frequencies
        cszⱼ = sqrt(gamma_z[s] * Psz[s] / ρⱼ)
        cspⱼ = sqrt(gamma_p[s] * Psp[s] / ρⱼ)

        # k dot v_drift
        kvdⱼ = kz * vdsz[s]  # kx*vdsx + kz*vdsz

        # dn ~ n & v (continuity equation)
        M[ind + 1, ind + 1] = kvdⱼ
        M[ind + 1, ind + 2] = kx * nⱼ
        M[ind + 1, ind + 4] = kz * nⱼ

        # dv ~ n & v (momentum equation)
        M[ind + 2, ind + 1] = kx * cspⱼ^2 / nⱼ
        M[ind + 4, ind + 1] = kz * cszⱼ^2 / nⱼ
        M[ind + 2, ind + 2] = kvdⱼ
        M[ind + 3, ind + 3] = kvdⱼ
        M[ind + 4, ind + 4] = kvdⱼ
        M[ind + 3, ind + 2] = -1im * wcⱼ
        M[ind + 2, ind + 3] = 1im * wcⱼ

        # dv ~ E (Lorentz force)
        M[ind + 2, SJ + 1] = im * qⱼ / mⱼ
        M[ind + 3, SJ + 2] = im * qⱼ / mⱼ
        M[ind + 4, SJ + 3] = im * qⱼ / mⱼ

        # dv ~ B (magnetic force and pressure anisotropy)
        ∆ⱼ = B0 == 0 ? 0 : (Psp[s] - Psz[s]) / B0
        M[ind + 2, SJ + 4] = kz * ∆ⱼ / ρⱼ
        M[ind + 2, SJ + 5] = -im * qⱼ / mⱼ * vdsz[s]
        M[ind + 4, SJ + 4] = kx * ∆ⱼ / ρⱼ
        M[ind + 3, SJ + 4] = 1im * qⱼ / mⱼ * vdsz[s]
        M[ind + 3, SJ + 5] = kz * ∆ⱼ / ρⱼ

        # dE ~ n
        M[SJ + 3, ind + 1] = -im * qⱼ * vdsz[s] / ε0
        # dE ~ v
        M[SJ + 1, ind + 2] = -im * qⱼ * nⱼ / ε0
        M[SJ + 2, ind + 3] = -im * qⱼ * nⱼ / ε0
        M[SJ + 3, ind + 4] = -im * qⱼ * nⱼ / ε0
    end

    # E(B): Faraday's law contribution to E equation
    c2 = c0 * c0
    M[SJ + 1, SJ + 5] = c2 * kz
    M[SJ + 2, SJ + 4] = -c2 * kz
    M[SJ + 2, SJ + 6] = c2 * kx
    M[SJ + 3, SJ + 5] = -c2 * kx

    # B(E): Faraday's law
    M[SJ + 4, SJ + 2] = -kz
    M[SJ + 5, SJ + 1] = kz
    M[SJ + 5, SJ + 3] = -kx
    M[SJ + 6, SJ + 2] = kx

    # Solve eigenvalue problem
    return eigen(M).values
end


function solve_fluid_dispersion(species, B0, kx, kz)
    params = create_fluid_params(species, B0)
    return solve_fluid_dispersion(params, kx, kz), params
end
