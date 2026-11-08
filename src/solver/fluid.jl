"""
    BOFluid(; gamma_perp = 1.0, gamma_para = 1.0)

Multi-fluid electromagnetic solver.

One polytropic index for every population, or a collection with one entry per population.
"""
@kwdef struct BOFluid{P, Z} <: AbstractDispersionAlgorithm
    gamma_perp::P = 1.0
    gamma_para::Z = 1.0
end

_per_species(x::Number, i) = x
_per_species(x, i) = x[i]

function _fluid_param(s::Population, gamma_p, gamma_z, ::Type{T}) where {T}
    vth_perp, vth_para = T.(vth_perp_para(s.vdf))
    q, m, n = T(s.particle.q), T(s.particle.m), T(s.n)
    return (;
        q, m, n, vdz = T(mean(s.vdf)[3]), Pz = n * m * vth_para^2 / 2, Pp = n * m * vth_perp^2 / 2,
        gamma_z = T(gamma_z), gamma_p = T(gamma_p),
    )
end

function prepare(alg::BOFluid, species, B0, ::Type{T}) where {T}
    S = length(species)
    for γ in (alg.gamma_perp, alg.gamma_para)
        γ isa Number || length(γ) == S ||
            throw(ArgumentError("BOFluid polytropic index $γ: need a scalar or one entry per population ($S)"))
    end
    return ntuple(i -> _fluid_param(species[i], _per_species(alg.gamma_perp, i), _per_species(alg.gamma_para, i), T), S)
end

# Matrix dimensions: 4 variables per species (n, vx, vy, vz) + 6 fields (Ex, Ey, Ez, Bx, By, Bz)
matrix_size(::BOFluid, species) = 4 * length(species) + 6

function _assemble_fluid_species!(M, s, ind, SJ, kx, kz, B0)
    q = s.q
    m = s.m
    n = s.n
    vdz = s.vdz
    Pz = s.Pz
    Pp = s.Pp

    ρ = m * n # Mass density (kg/m³)
    wc = q * B0 / m # Cyclotron frequency
    csz = sqrt(s.gamma_z * Pz / ρ)
    csp = sqrt(s.gamma_p * Pp / ρ)

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

function dispersion_matrix!(M, pb::DispersionProblem, alg::BOFluid; c2 = c0^2)
    B0, kx, kz = pb.B0, pb.kx, pb.kz
    species = pb.species
    S = length(species)
    SJ = 4 * S
    for i in 1:S
        _assemble_fluid_species!(M, species[i], (i - 1) * 4, SJ, kx, kz, B0)
    end
    _B_E_part!(M, SJ, kx, kz; c2)
    return M
end

