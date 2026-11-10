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

# Pressure anisotropy pulls along the perturbed field line, a force ∝ ΔP δB/B₀: no B₀-free limit.
function _fluid_param(s, gamma_p, gamma_z, ::Type{T}) where {T}
    c = T(s.c)
    vth_perp, vth_para = T.(vth_perp_para(s.vdf)) ./ c
    wc, wp = T(s.Omega), sqrt(abs(T(s.Pi2)))
    Δ = (vth_para^2 - vth_perp^2) / 2
    aniso = iszero(Δ) || iszero(wp) ? zero(T) : iszero(wc) ?
        throw(ArgumentError("BOFluid: an anisotropic population needs B0 ≠ 0")) : Δ * wp / wc
    return (;
        wc, wp, sgn = sign(T(s.Pi2)), vdz = T(mean(s.vdf)[3]) / c, aniso,
        csp2 = T(gamma_p) * vth_perp^2 / 2, csz2 = T(gamma_z) * vth_para^2 / 2,
    )
end

function prepare(alg::BOFluid, populations, ::Type{T}) where {T}
    S = length(populations)
    for γ in (alg.gamma_perp, alg.gamma_para)
        γ isa Number || length(γ) == S ||
            throw(ArgumentError("BOFluid polytropic index $γ: need a scalar or one entry per population ($S)"))
    end
    return ntuple(i -> _fluid_param(populations[i], _per_species(alg.gamma_perp, i), _per_species(alg.gamma_para, i), T), S)
end

# Matrix dimensions: 4 variables per species (n, vx, vy, vz) + 6 fields (Ex, Ey, Ez, Bx, By, Bz)
matrix_size(::BOFluid, species) = 4 * length(species) + 6

# Fields are E and B with c = 1; a species' δn/n and v (in c) are scaled by ω_p/Ω, so that its
# couplings to the fields need ω_p and Ω only, never their ratio, and survive B₀ = 0.
function _assemble_fluid_species!(M, s, ind, SJ, kx, kz)
    (; wc, wp, vdz, aniso, csp2, csz2) = s
    kvd = kz * vdz

    # dn ~ v (continuity equation)
    M[ind + 1, ind + 1] = kvd
    M[ind + 1, ind + 2] = kx
    M[ind + 1, ind + 4] = kz

    # dv ~ n & v (momentum equation)
    M[ind + 2, ind + 1] = kx * csp2
    M[ind + 4, ind + 1] = kz * csz2
    M[ind + 2, ind + 2] = kvd
    M[ind + 3, ind + 3] = kvd
    M[ind + 4, ind + 4] = kvd
    M[ind + 3, ind + 2] = -1im * wc
    M[ind + 2, ind + 3] = 1im * wc

    # dv ~ E (Lorentz force)
    M[ind + 2, SJ + 1] = 1im * wp
    M[ind + 3, SJ + 2] = 1im * wp
    M[ind + 4, SJ + 3] = 1im * wp

    # dv ~ B (magnetic force and pressure anisotropy)
    M[ind + 2, SJ + 4] = kz * aniso
    M[ind + 2, SJ + 5] = -1im * wp * vdz
    M[ind + 4, SJ + 4] = kx * aniso
    M[ind + 3, SJ + 4] = 1im * wp * vdz
    M[ind + 3, SJ + 5] = kz * aniso

    # dE ~ n & v (current density)
    a = s.sgn * wp
    M[SJ + 3, ind + 1] = -1im * a * vdz
    M[SJ + 1, ind + 2] = -1im * a
    M[SJ + 2, ind + 3] = -1im * a
    M[SJ + 3, ind + 4] = -1im * a
    return nothing
end

function dispersion_matrix!(M, species, k::Wavevector, alg::BOFluid)
    kx, kz = k.k_perp, k.k_para
    S = length(species)
    SJ = 4 * S
    for i in 1:S
        _assemble_fluid_species!(M, species[i], (i - 1) * 4, SJ, kx, kz)
    end
    _B_E_part!(M, SJ, kx, kz)
    return M
end
