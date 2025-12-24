# Electromagnetic fluid dispersion relation solver

"""
    FluidSpecies{T}

Fluid species parameters for the multi-fluid dispersion relation solver.

# Fields
- `q`: Charge
- `m`: Mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity (in units of c)
- `gamma_z`: Parallel polytrope exponent (default: 1.0 for isothermal)
- `gamma_p`: Perpendicular polytrope exponent (default: 1.0 for isothermal)
"""
struct FluidSpecies{T <: Number}
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    gamma_z::T     # Parallel polytrope exponent
    gamma_p::T     # Perpendicular polytrope exponent
    q::T           # Charge
    m::T           # Mass
end

"""
    FluidSpecies(n, Tz, Tp = Tz; vdz=0.0, gamma_z=1.0, gamma_p=1.0, particle=:p, ...)
"""
function FluidSpecies(n, Tz = 0.0, Tp = Tz; vdz = 0.0, gamma_z = 1.0, gamma_p = 1.0, Z = nothing, A = nothing, q = nothing, m = nothing, particle = :p)
    q, m = _charge_mass(particle, Z, A, q, m)
    return FluidSpecies(promote(n, Tz, Tp, vdz, gamma_z, gamma_p, q, m)...)
end

FluidSpecies(p::ParticleLike, args...; kw...) = FluidSpecies(args...; kw...)

plasma_frequency(s::FluidSpecies) = plasma_frequency(s.q, s.n, s.m)

# Matrix dimensions: 4 variables per species (n, vx, vy, vz) + 6 fields (Ex, Ey, Ez, Bx, By, Bz)
_fluid_size(S::Int) = 4 * S + 6
_fluid_size(species) = _fluid_size(length(species))

_gamma_z(s::FluidSpecies) = s.gamma_z
_gamma_p(s::FluidSpecies) = s.gamma_p
_gamma_z(s) = 1.0
_gamma_p(s) = 1.0
# Pressures: P = n * k_B * T, with T in Kelvin (Tz in eV → T = Tz * q / k_B)
_Pz(s) = s.n * s.Tz * qe
_Pp(s) = s.n * s.Tp * qe

function _assemble_fluid_species!(M, s, ind, SJ, kx, kz, B0)
    q = s.q
    m = s.m
    n = s.n
    vdz = s.vdz * c0
    Pz = _Pz(s)
    Pp = _Pp(s)

    ρ = m * n # Mass density (kg/m³)
    wc = q * B0 / m # Cyclotron frequency
    csz = sqrt(_gamma_z(s) * Pz / ρ)
    csp = sqrt(_gamma_p(s) * Pp / ρ)

    # k dot v_drift
    kvd = kz * vdz

    # dn ~ n & v (continuity equation)
    M[ind + 1, ind + 1] = kvd
    M[ind + 1, ind + 2] = kx * n
    M[ind + 1, ind + 4] = kz * n

    # dv ~ n & v (momentum equation)
    M[ind + 2, ind + 1] = kx * csp^2 / n
    M[ind + 4, ind + 1] = kz * csz^2 / n
    M[ind + 2, ind + 2] = kvd
    M[ind + 3, ind + 3] = kvd
    M[ind + 4, ind + 4] = kvd
    M[ind + 3, ind + 2] = -1im * wc
    M[ind + 2, ind + 3] = 1im * wc

    # dv ~ E (Lorentz force)
    qm = q / m
    M[ind + 2, SJ + 1] = 1im * qm
    M[ind + 3, SJ + 2] = 1im * qm
    M[ind + 4, SJ + 3] = 1im * qm

    # dv ~ B (magnetic force and pressure anisotropy)
    Δ = B0 == 0 ? 0.0 : (Pp - Pz) / B0 / ρ
    M[ind + 2, SJ + 4] = kz * Δ
    M[ind + 2, SJ + 5] = -1im * qm * vdz
    M[ind + 4, SJ + 4] = kx * Δ
    M[ind + 3, SJ + 4] = 1im * qm * vdz
    M[ind + 3, SJ + 5] = kz * Δ

    # dE ~ n
    M[SJ + 3, ind + 1] = -1im * q * vdz / ε0
    # dE ~ v
    M[SJ + 1, ind + 2] = -1im * q * n / ε0
    M[SJ + 2, ind + 3] = -1im * q * n / ε0
    M[SJ + 3, ind + 4] = -1im * q * n / ε0
    return nothing
end

function build_fluid_dispersion_matrix(species, args...; kw...)
    NN = _fluid_size(species)
    M = zeros(ComplexF64, NN, NN)
    return build_fluid_dispersion_matrix!(M, species, args...; kw...)
end

# Build the fluid dispersion matrix in-place.
function build_fluid_dispersion_matrix!(M, species, kx, kz, B0; c2 = c0^2)
    S = length(species)
    SJ = 4 * S
    for i in 1:S
        ind = (i - 1) * 4
        _assemble_fluid_species!(M, species[i], ind, SJ, kx, kz, B0)
    end
    _B_E_part!(M, SJ, kx, kz; c2)
    return M
end

"""
    solve_fluid_dispersion(species, B0, kx, kz)

Solve the multi-fluid electromagnetic dispersion relation using matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

The fluid model includes:
- Continuity equation: ∂n/∂t + ∇·(nv) = 0
- Momentum equation: m(∂v/∂t + v·∇v) = q(E + v×B) - ∇P/n
- Maxwell's equations for E and B

See also: [`FluidSpecies`](@ref)
"""
function solve_fluid_dispersion(species, B0, kx, kz)
    M = build_fluid_dispersion_matrix(species, kx, kz, B0)
    return eigvals!(M)
end

function solve_fluid_dispersion!(M, species, kx, kz, B0)
    fill!(M, zero(eltype(M)))
    build_fluid_dispersion_matrix!(M, species, kx, kz, B0)
    return eigvals!(M)
end

function solve_fluid_dispersion(species, B0, ks::AbstractVector, θ)
    NN = _fluid_size(species)
    M = zeros(ComplexF64, NN, NN)
    ωs = @showprogress desc = "Solving fluid dispersion..." map(ks) do k
        kx, kz = k .* sincos(θ)
        solve_fluid_dispersion!(M, species, kx, kz, B0)
    end
    return DispersionSolution(ks, ωs)
end
