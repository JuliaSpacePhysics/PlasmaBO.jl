"""
    BOPBK(; N = 2)

Dispersion solver using the PBK matrix formulation.

`N` controls the truncation order of the cyclotron harmonic index used to build
the dispersion matrix.
"""
@kwdef struct BOPBK <: AbstractDispersionAlgorithm
    N::Int = 2
end

# Past this κ∥, roundoff spreads the pole's chain into the roots (docs/src/math_notes.md).
const _PBK_KAPPA_PARA_MAX = 24

function PBK_param(s, ::Type{T}) where {T}
    iszero(s.Omega) && throw(ArgumentError("BOPBK expands about the gyration: needs B0 ≠ 0"))
    d = s.vdf
    d isa ProductBiKappa || throw(ArgumentError("BOPBK solves ProductBiKappa species, got $(nameof(typeof(d)))"))
    isinteger(d.kappa_para) && 1 ≤ d.kappa_para ≤ _PBK_KAPPA_PARA_MAX || throw(ArgumentError(
        "BOPBK needs an integer κ∥ in 1:$_PBK_KAPPA_PARA_MAX, beyond which roundoff moves its roots, got $(d.kappa_para)"
    ))
    c = T(s.c)
    vtx, vtz = T.(kappa_widths(d)) ./ c
    wc = T(s.Omega)
    return (;
        u0 = T(d.vd_para) / c, sigma = T(d.loss_cone), κz = Int(d.kappa_para), κx = T(d.kappa_perp),
        wp2 = T(s.Pi2), wc, vtz, vtx, rhoc = vtx / abs(wc),   # signed wp2: the response is linear in n
    )
end

prepare(::BOPBK, populations, ::Type{T}) where {T} = map(sp -> PBK_param(sp, T), populations)

# Pole csn = n ωc + k∥ ξ
_ξ(p) = p.u0 - 1im * sqrt(p.κz) * p.vtz

# Ez couples through tanθ = k⊥/k∥. Every such term carries enough powers of k∥ to cancel it, so
# the coefficients below take `r` factors of k∥ already divided out, which keeps k∥ = 0 regular.
function _csl(p, l, kz, r = 0)
    κz = p.κz
    z = 2im * sqrt(κz) * p.vtz
    return exp(loggamma(κz + 1) + loggamma(2 * κz - l + 2) - loggamma(κz - l + 2) - loggamma(2 * κz + 1)) *
        z^r * (z * kz)^(l - 1 - r)
end

bsnl(p, n, l, kz, r) = -n * p.wc * _csl(p, l, kz, r)
bsl(p, l, kz, r) = -l * kz^(2 - r) * p.vtx^2 * _csl(p, l, kz) / 2

# J_ν(z) / (z/2)^ν, finite at z = 0
function _besselj_scaled(ν, z::T) where {T}
    t = (z / 2)^2
    t < sqrt(eps(T)) && return (1 - t / (ν + 1)) / gamma(T(ν + 1))
    return besselj(ν, z) / (z / 2)^ν
end

# J_ν'(z) / (z/2)^|ν - 1|
function _dbesselj_scaled(ν, z)
    iszero(ν) && return -_besselj_scaled(1, z)
    return (_besselj_scaled(ν - 1, z) - (z / 2)^2 * _besselj_scaled(ν + 1, z)) / 2
end

# (p1, p2, p3, p4, m) of the S integrals, in `_pbk_S`'s order; only the first three enter at ν = 0.
const _PBK_S_TERMS = (
    (2, 0, 1, 1, 2), (1, 1, 2, 1, 0), (0, 2, 3, 1, 0), (2, 0, 1, 1, 0), (2, 0, 1, 2, 0),
    (2, 0, -1, 1, 2), (1, 1, 2, 2, 0), (1, 1, 0, 1, 2), (0, 2, 3, 2, 0), (0, 2, 1, 1, 2),
)

# aᵐ S per term, a = k⊥ρc. As a → 0, S ~ a^(p3 - 3 + p1ν + p2|ν - 1|), ν = |n|: some S diverge, but
# each such S is multiplied by λ or k⊥² (folded in through `m`) or by a factor vanishing at n = 0.
# The integrand carries the Bessel powers factored out, so the folded result is finite at k⊥ = 0.
# p1 + p2 = 2 makes S even in n. One quadrature for all terms shares the Bessel evaluations.
function funcS_pbk(p, ν, λ::T) where {T}
    κx = p.κx
    sgm = p.sigma
    eps0 = eps(one(T))
    C = 4 * κx^(-sgm - 1) * exp(loggamma(κx + sgm + 1) - loggamma(κx)) / gamma(sgm + 1)
    nterms = iszero(ν) ? 3 : length(_PBK_S_TERMS)
    a = sqrt(2 * λ)
    eJ(p1, p2) = p1 * ν + p2 * abs(ν - 1)
    I = quadgk(zero(T), T(Inf)) do x
        jn, djn = _besselj_scaled(ν, x * a), _dbesselj_scaled(ν, x * a)
        SVector(ntuple(length(_PBK_S_TERMS)) do i
            p1, p2, p3, p4 = _PBK_S_TERMS[i]
            i > nterms && return zero(T)
            num = (x + eps0)^(2 * sgm + p3) * x^eJ(p1, p2) * jn^p1 * djn^p2
            # log1p: the denominator itself overflows Double64 (max ~1e308)
            num * exp(-(κx + sgm + p4) * log1p(x^2 / κx))
        end)
    end[1]
    return ntuple(length(_PBK_S_TERMS)) do i
        p1, p2, p3, _, m = _PBK_S_TERMS[i]
        C * a^(p3 - 3 + m + eJ(p1, p2)) * I[i] / 2^eJ(p1, p2)
    end
end

# Combinations of the S integrals for harmonics ±ν, shared by every l. At ν = 0 those entering only
# with a factor n are left zero: the products vanish, while some of the integrals diverge as k⊥ → 0.
function _pbk_S(p, ν, kx)
    λ = kx^2 * p.rhoc^2 / 2
    c, σ = (p.κx + p.sigma + 1) / p.κx, p.sigma
    kx2s2, s4, s6, s2, s17a, s17b, s38a, s38b, s59a, s59b = funcS_pbk(p, ν, λ)
    z = zero(λ)
    iszero(ν) && return (; s17 = z, s2 = z, kx2s2 = kx2s2 / p.rhoc^2, s38 = z, s4, s59 = z, s6)
    return (;
        s17 = c * s17a - σ * s17b, s2, kx2s2 = kx2s2 / p.rhoc^2,
        s38 = c * s38a - σ * s38b, s4, s59 = c * s59a - σ * s59b, s6,
    )
end

b1snl(S, p, n, l, kz; r = 0) = S.s17 * bsnl(p, n, l, kz, r)
b2snl(S, p, l, kz; r = 0) = S.s2 * bsl(p, l, kz, r)
# k⊥² b2, finite at k⊥ = 0 where b2 at n = 0 is not
kx2b2snl(S, p, l, kz; r = 0) = S.kx2s2 * bsl(p, l, kz, r)
b3snl(S, p, n, l, kz; r = 0) = S.s38 * bsnl(p, n, l, kz, r)
b4snl(S, p, l, kz; r = 0) = S.s4 * bsl(p, l, kz, r)
b5snl(S, p, n, l, kz) = S.s59 * bsnl(p, n, l, kz, 0)
b6snl(S, p, l, kz) = S.s6 * bsl(p, l, kz, 0)

bx11snl(wp2, n, b1) = -1im * wp2 * n^2 * b1
bx12snl(wp2, n, b2) = -1im * wp2 * n^2 * b2
bx21snl(wp2, n, b3) = wp2 * n * b3
bx22snl(wp2, n, b4) = wp2 * n * b4
bx31snl(wp2, wc, n, kx, ξ, b1, b2_kz) = -1im * kx * wp2 * n * (b1 * ξ + b2_kz) / wc
bx32snl(wp2, wc, n, kx, ξ, b2) = -1im * kx * wp2 * n * b2 * ξ / wc
bx33snl(wp2, wc, n, kx, b1_kz) = -1im * kx * wp2 * n * b1_kz / wc

by21snl(wp2, b5) = -1im * wp2 * b5
by22snl(wp2, b6) = -1im * wp2 * b6
by31snl(wp2, wc, kx, ξ, b3, b4_kz) = -kx * wp2 * (b3 * ξ + b4_kz) / wc
by32snl(wp2, wc, kx, ξ, b4) = -kx * wp2 * b4 * ξ / wc
by33snl(wp2, wc, kx, b3_kz) = -kx * wp2 * b3_kz / wc

bz31snl(wp2, wc, kx, ξ, b1, kx2b2_kz) = -1im * wp2 * (kx^2 * b1 * ξ^2 + 2 * kx2b2_kz * ξ) / wc^2
bz32snl(wp2, wc, ξ, kx2b2) = -1im * wp2 * kx2b2 * ξ^2 / wc^2
bz33snl(wp2, wc, kx, ξ, b1_kz, kx2b2_kz2) = -1im * wp2 * (2 * kx^2 * b1_kz * ξ + kx2b2_kz2) / wc^2
bz34snl(wp2, wc, kx, b1_kz2) = -1im * kx^2 * wp2 * b1_kz2 / wc^2

# The E-couplings of level l of harmonic n, as `_coupling!` takes them: the head, at (ω - csn)^-l,
function _pbk_head(S, p, n, l, kx, kz, ξ)
    (; κz, wp2, wc) = p
    b1, b3 = b1snl(S, p, n, l, kz), b3snl(S, p, n, l, kz)
    p13 = bx31snl(wp2, wc, n, kx, ξ, b1, b2snl(S, p, l, kz; r = 1))
    p23 = by31snl(wp2, wc, kx, ξ, b3, b4snl(S, p, l, kz; r = 1))
    p33 = bz31snl(wp2, wc, kx, ξ, b1, kx2b2snl(S, p, l, kz; r = 1))
    if l <= κz
        b1_kz = b1snl(S, p, n, l + 1, kz; r = 1)
        p13 += bx33snl(wp2, wc, n, kx, b1_kz)
        p23 += by33snl(wp2, wc, kx, b3snl(S, p, n, l + 1, kz; r = 1))
        p33 += bz33snl(wp2, wc, kx, ξ, b1_kz, kx2b2snl(S, p, l + 1, kz; r = 2))
    end
    l <= κz - 1 && (p33 += bz34snl(wp2, wc, kx, b1snl(S, p, n, l + 2, kz; r = 2)))
    return (bx11snl(wp2, n, b1), bx21snl(wp2, n, b3), p13, by21snl(wp2, b5snl(S, p, n, l, kz)), p23, p33)
end

# and the tail, at (ω - csn)^-(l+1)
function _pbk_tail(S, p, n, l, kx, kz, ξ)
    (; wp2, wc) = p
    b2, b4 = b2snl(S, p, l, kz), b4snl(S, p, l, kz)
    return (
        bx12snl(wp2, n, b2), bx22snl(wp2, n, b4), bx32snl(wp2, wc, n, kx, ξ, b2),
        by22snl(wp2, b6snl(S, p, l, kz)), by32snl(wp2, wc, kx, ξ, b4), bz32snl(wp2, wc, ξ, kx2b2snl(S, p, l, kz)),
    )
end

# One chain of κ∥ + 2 states per harmonic, scaled by the pole width s; see docs/src/math_notes.md.
# Each component (x, y, z) has the chains of every species and its current J last.
_pbk_nstates(p, N) = (2N + 1) * (p.κz + 2)

matrix_size(alg::BOPBK, params) = 3 * (sum(p -> _pbk_nstates(p, alg.N), params) + 1) + 6

# Fills the chains of one species after row `r` of each component block of `L` rows; returns its last row.
function _pbk_species!(M, p, r, L, kx, kz, N)
    ξ = _ξ(p)
    s = 2 * sqrt(p.κz) * p.vtz * abs(kz)
    iszero(s) && (s = one(s))
    m = p.κz + 2
    Ss = [_pbk_S(p, ν, kx) for ν in 0:N]
    for n in -N:N
        S = Ss[abs(n) + 1]
        csn = n * p.wc + kz * ξ
        for c in 0:2
            M[(c + 1) * L, c * L + r + 1] = 1
        end
        for j in 1:m
            r += 1
            for c in 0:2
                M[c * L + r, c * L + r] = csn
                j < m && (M[c * L + r, c * L + r + 1] = s)
            end
            b = j < m ? _pbk_head(S, p, n, j, kx, kz, ξ) : ntuple(_ -> zero(csn), 6)
            j > 1 && (b = b .+ _pbk_tail(S, p, n, j - 1, kx, kz, ξ))
            _coupling!(view(M, r:L:(r + 2L), (3L + 1):(3L + 3)), b ./ s^(j - 1))
        end
    end
    return r
end

# The closed form of Z_κ∥ follows the Landau contour only for k∥ > 0. k∥ < 0 is the mirror image
# z → -z with drifts reversed, mapped back by flipping the signs of Jz's chains, Ez, Bx and By.
function dispersion_matrix!(M, params, k::Wavevector, alg::BOPBK)
    T = real(eltype(M))
    kx, kz = T(k.k_perp), T(k.k_para)
    kz < 0 || (_pbk_matrix!(M, params, kx, kz, alg.N); return M)
    L = _pbk_matrix!(M, map(p -> merge(p, (; u0 = -p.u0)), params), kx, -kz, alg.N)
    for I in ((2L + 1):3L, (3L + 3):(3L + 5))
        M[I, :] .*= -1
        M[:, I] .*= -1
    end
    return M
end

function _pbk_matrix!(M, params, kx, kz, N)
    L = sum(p -> _pbk_nstates(p, N), params) + 1
    r = 0
    for p in params
        r = _pbk_species!(M, p, r, L, kx, kz, N)
    end
    iwp2 = 1im * sum(p -> p.wp2, params)
    for c in 1:3
        M[c * L, 3L + c] = iwp2
        M[3L + c, c * L] = -1im
    end
    _B_E_part!(M, 3L, kx, kz)
    return L
end
