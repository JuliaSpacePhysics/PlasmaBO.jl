module PBK
using QuadGK: quadgk
using ..PlasmaBO: AbstractDispersionAlgorithm, DispersionProblem, _realtype, c0, ε0
using ..PlasmaBO: besselj, loggamma, gamma
using PlasmaPopulations: Population, ProductBiKappa, kappa_widths
import ..PlasmaBO: dispersion_matrix!, matrix_size, prepare

"""
    BOPBK(; N = 2)

Dispersion solver using the PBK matrix formulation.

`N` controls the truncation order of the cyclotron harmonic index used to build
the dispersion matrix.
"""
@kwdef struct BOPBK <: AbstractDispersionAlgorithm
    N::Int = 2
end

function PBK_param(s::Population, B0, ::Type{T}) where {T}
    d = s.vdf
    d isa ProductBiKappa || throw(ArgumentError("BOPBK solves ProductBiKappa species, got $(nameof(typeof(d)))"))
    isinteger(d.kappa_para) || throw(ArgumentError("BOPBK needs an integer κ∥, got $(d.kappa_para)"))
    m, q = T(s.particle.m), T(s.particle.q)
    wp2 = T(s.n) * q^2 / (m * T(ε0))   # signed: the response is linear in n
    wc = (T(B0) * q) / m
    vtx, vtz = T.(kappa_widths(d))
    return (;
        u0 = T(d.vd_para), sigma = T(d.loss_cone), κz = Int(d.kappa_para), κx = T(d.kappa_perp),
        wp2, wc, vtz, vtx, rhoc = vtx / abs(wc),
    )
end

prepare(::BOPBK, species, B0, ::Type{T}) where {T} = map(sp -> PBK_param(sp, B0, T), species)

_pbk_per_n(κ) = (κ + 1) * (κ + 4) ÷ 2

function _pbk_len_total(species, Ns)
    len = 0
    for sp in species
        len += (2 * Ns + 1) * _pbk_per_n(Int(sp.κz))
    end
    return len
end

_pbk_len_sub(species, N) = _pbk_len_total(species, N) + 1

function _pbk_species_offset(species, s, Ns)
    off = 0
    for i in 1:(s - 1)
        sp = species[i]
        off += (2 * Ns + 1) * _pbk_per_n(Int(sp.κz))
    end
    return off
end

function _pbk_snl_index(species, s, idx_n, l, j, Ns)
    sp = species[s]
    idx_n < 1 || idx_n > 2 * Ns + 1 && error("idx_n out of range")
    l < 1 || l > Int(sp.κz) + 1 && error("l out of range")
    j < 1 || j > l + 1 && error("j out of range")

    per_n = _pbk_per_n(Int(sp.κz))
    off = _pbk_species_offset(species, s, Ns)
    off_n = (idx_n - 1) * per_n
    off_l = (l + 2) * (l - 1) ÷ 2
    return off + off_n + off_l + j
end

function _pbk_block_first_index(len_sub, MatrixNo)
    MatrixNo == 1 && return 0
    MatrixNo == 2 && return len_sub
    MatrixNo == 3 && return 2 * len_sub
    error("MatrixNo must be 1, 2, or 3")
end

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

function funcS_pbk(p, n, λ::T, p1, p2, p3, p4; EPS0 = 1.0e-2) where {T}
    κx = p.κx
    sgm = p.sigma
    eps0 = eps(one(T))

    # Denominator is computed via `log1p` to avoid exponent overflow in
    # extended precision (Double64 maxes at ~1e308 like Float64).
    if λ > EPS0
        S0 = 4 * (2 * λ)^(-sgm - 2) * κx^(-sgm - 1) * exp(loggamma(κx + sgm + 1) - loggamma(κx)) / gamma(sgm + 1)
        function integrand1(x)
            jn = besselj(n, x)
            djn = (besselj(n - 1, x) - besselj(n + 1, x)) / 2
            num = (x + eps0)^(2 * sgm + p3) * jn^p1 * djn^p2
            log_den = (κx + sgm + p4) * log1p(x^2 / λ / κx / 2)
            return num * exp(-log_den)
        end

        tmp = quadgk(integrand1, zero(T), T(Inf))[1]
        return S0 * tmp
    else
        S0 = 4 * (2 * λ)^(p3 / 2 - 3//2) * κx^(-sgm - 1) * exp(loggamma(κx + sgm + 1) - loggamma(κx)) / gamma(sgm + 1)
        a = sqrt(2 * λ)
        function integrand2(x)
            jn = besselj(n, x * a)
            djn = (besselj(n - 1, x * a) - besselj(n + 1, x * a)) / 2
            num = (x + eps0)^(2 * sgm + p3) * jn^p1 * djn^p2
            log_den = (κx + sgm + p4) * log1p(x^2 / κx)
            return num * exp(-log_den)
        end
        tmp = quadgk(integrand2, zero(T), T(Inf))[1]
        return S0 * tmp
    end
end

function b1snl(p, n, l, kz, λ; r = 0, EPS0 = 1.0e-2)
    S1 = funcS_pbk(p, n, λ, 2, 0, 1, 2; EPS0)
    S7 = funcS_pbk(p, n, λ, 2, 0, -1, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S17 = (κx + sgm + 1) / κx * S1 - 2 * sgm * λ * S7
    return S17 * bsnl(p, n, l, kz, r)
end

function b2snl(p, n, l, kz, λ; r = 0, EPS0 = 1.0e-2)
    S2 = funcS_pbk(p, n, λ, 2, 0, 1, 1; EPS0)
    return S2 * bsl(p, l, kz, r)
end

function b3snl(p, n, l, kz, λ; r = 0, EPS0 = 1.0e-2)
    S3 = funcS_pbk(p, n, λ, 1, 1, 2, 2; EPS0)
    S8 = funcS_pbk(p, n, λ, 1, 1, 0, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S38 = (κx + sgm + 1) / κx * S3 - 2 * sgm * λ * S8
    return S38 * bsnl(p, n, l, kz, r)
end

function b4snl(p, n, l, kz, λ; r = 0, EPS0 = 1.0e-2)
    S4 = funcS_pbk(p, n, λ, 1, 1, 2, 1; EPS0)
    return S4 * bsl(p, l, kz, r)
end

function b5snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S5 = funcS_pbk(p, n, λ, 0, 2, 3, 2; EPS0)
    S9 = funcS_pbk(p, n, λ, 0, 2, 1, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S59 = (κx + sgm + 1) / κx * S5 - 2 * sgm * λ * S9
    return S59 * bsnl(p, n, l, kz, 0)
end

function b6snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S6 = funcS_pbk(p, n, λ, 0, 2, 3, 1; EPS0)
    return S6 * bsl(p, l, kz, 0)
end

bx11snl(wp2, n, b1) = -1im * ε0 * wp2 * n^2 * b1
bx12snl(wp2, n, b2) = -1im * ε0 * wp2 * n^2 * b2
bx21snl(wp2, n, b3) = ε0 * wp2 * n * b3
bx22snl(wp2, n, b4) = ε0 * wp2 * n * b4
bx31snl(wp2, wc, n, kx, ξ, b1, b2_kz) = -1im * ε0 * kx * wp2 * n * (b1 * ξ + b2_kz) / wc
bx32snl(wp2, wc, n, kx, ξ, b2) = -1im * ε0 * kx * wp2 * n * b2 * ξ / wc
bx33snl(wp2, wc, n, kx, b1_kz) = -1im * ε0 * kx * wp2 * n * b1_kz / wc

by21snl(wp2, b5) = -1im * ε0 * wp2 * b5
by22snl(wp2, b6) = -1im * ε0 * wp2 * b6
by31snl(wp2, wc, kx, ξ, b3, b4_kz) = -ε0 * kx * wp2 * (b3 * ξ + b4_kz) / wc
by32snl(wp2, wc, kx, ξ, b4) = -ε0 * kx * wp2 * b4 * ξ / wc
by33snl(wp2, wc, kx, b3_kz) = -ε0 * kx * wp2 * b3_kz / wc

bz31snl(wp2, wc, kx, ξ, b1, b2_kz) = -1im * ε0 * kx^2 * wp2 * (b1 * ξ^2 + 2 * b2_kz * ξ) / wc^2
bz32snl(wp2, wc, kx, ξ, b2) = -1im * ε0 * kx^2 * wp2 * b2 * ξ^2 / wc^2
bz33snl(wp2, wc, kx, ξ, b1_kz, b2_kz2) = -1im * ε0 * kx^2 * wp2 * (2 * b1_kz * ξ + b2_kz2) / wc^2
bz34snl(wp2, wc, kx, b1_kz2) = -1im * ε0 * kx^2 * wp2 * b1_kz2 / wc^2

lmax(κz; threshold = 10) = min(κz + 1, threshold)

function _pbk_add_species_Mxy!(
        M, params, s,
        kx, kz;
        MatrixNo, ExNo, EyNo, EzNo, ExyNo,
        bx10_by20, EPS0, len_sub, N
    )
    p = params[s]
    firstIndex = _pbk_block_first_index(len_sub, MatrixNo)
    λ = kx^2 * (p.rhoc^2) / 2
    κz = p.κz

    wp2 = p.wp2
    wc = p.wc
    ξ = _ξ(p)

    for (idx_n, n) in enumerate(-N:N)
        csn = n * wc + kz * ξ
        for l in 1:lmax(κz), jj in 1:(l + 1)
            snlj = _pbk_snl_index(params, s, idx_n, l, jj, N)

            M[snlj, firstIndex + snlj] += csn
            if jj < l + 1
                M[snlj, firstIndex + snlj + 1] += 1
            end

            if jj == l
                if MatrixNo == 1
                    b1 = b1snl(p, n, l, kz, λ; EPS0)
                    b2_kz = b2snl(p, n, l, kz, λ; r = 1, EPS0)
                    b3 = b3snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += bx11snl(wp2, n, b1)
                    M[snlj, end - EyNo] += bx21snl(wp2, n, b3)
                    M[snlj, end - EzNo] += bx31snl(wp2, wc, n, kx, ξ, b1, b2_kz)

                    if l <= κz
                        b1_kz = b1snl(p, n, l + 1, kz, λ; r = 1, EPS0)
                        M[snlj, end - EzNo] += bx33snl(wp2, wc, n, kx, b1_kz)
                    end
                else
                    b3 = b3snl(p, n, l, kz, λ; EPS0)
                    b4_kz = b4snl(p, n, l, kz, λ; r = 1, EPS0)
                    b5 = b5snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += -bx21snl(wp2, n, b3)
                    M[snlj, end - EyNo] += by21snl(wp2, b5)
                    M[snlj, end - EzNo] += by31snl(wp2, wc, kx, ξ, b3, b4_kz)

                    if l <= κz
                        b3_kz = b3snl(p, n, l + 1, kz, λ; r = 1, EPS0)
                        M[snlj, end - EzNo] += by33snl(wp2, wc, kx, b3_kz)
                    end
                end
            elseif jj == l + 1
                if MatrixNo == 1
                    b2 = b2snl(p, n, l, kz, λ; EPS0)
                    b4 = b4snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += bx12snl(wp2, n, b2)
                    M[snlj, end - EyNo] += bx22snl(wp2, n, b4)
                    M[snlj, end - EzNo] += bx32snl(wp2, wc, n, kx, ξ, b2)
                else
                    b4 = b4snl(p, n, l, kz, λ; EPS0)
                    b6 = b6snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += -bx22snl(wp2, n, b4)
                    M[snlj, end - EyNo] += by22snl(wp2, b6)
                    M[snlj, end - EzNo] += by32snl(wp2, wc, kx, ξ, b4)
                end
            end

        end
    end

    for idx_n in 1:(2 * N + 1), l in 1:(κz + 1)
        snl1 = _pbk_snl_index(params, s, idx_n, l, 1, N)
        M[len_sub, firstIndex + snl1] += 1
    end

    if s == length(params)
        M[len_sub, end - ExyNo] += bx10_by20
    end

    return nothing
end

function _pbk_add_species_Mz!(
        M, params, s, kx, kz;
        MatrixNo, ExNo, EyNo, EzNo,
        by20, EPS0, len_sub, N,
    )
    p = params[s]
    firstIndex = _pbk_block_first_index(len_sub, MatrixNo)
    λ = kx^2 * (p.rhoc^2) / 2
    κz = p.κz
    wp2 = p.wp2
    wc = p.wc
    ξ = _ξ(p)
    for (idx_n, n) in enumerate(-N:N)
        csn = n * wc + kz * ξ
        for l in 1:lmax(κz), jj in 1:(l + 1)
            snlj = _pbk_snl_index(params, s, idx_n, l, jj, N)

            M[snlj, firstIndex + snlj] += csn
            if jj < l + 1
                M[snlj, firstIndex + snlj + 1] += 1
            end

            if jj == l
                b1 = b1snl(p, n, l, kz, λ; EPS0)
                b2_kz = b2snl(p, n, l, kz, λ; r = 1, EPS0)
                b3 = b3snl(p, n, l, kz, λ; EPS0)
                b4_kz = b4snl(p, n, l, kz, λ; r = 1, EPS0)

                M[snlj, end - ExNo] += bx31snl(wp2, wc, n, kx, ξ, b1, b2_kz)
                M[snlj, end - EyNo] += -by31snl(wp2, wc, kx, ξ, b3, b4_kz)
                M[snlj, end - EzNo] += bz31snl(wp2, wc, kx, ξ, b1, b2_kz)

                if l <= κz
                    b1_kz = b1snl(p, n, l + 1, kz, λ; r = 1, EPS0)
                    b2_kz2 = b2snl(p, n, l + 1, kz, λ; r = 2, EPS0)
                    b3_kz = b3snl(p, n, l + 1, kz, λ; r = 1, EPS0)

                    M[snlj, end - ExNo] += bx33snl(wp2, wc, n, kx, b1_kz)
                    M[snlj, end - EyNo] += -by33snl(wp2, wc, kx, b3_kz)
                    M[snlj, end - EzNo] += bz33snl(wp2, wc, kx, ξ, b1_kz, b2_kz2)
                end

                if l <= κz - 1
                    b1_kz2 = b1snl(p, n, l + 2, kz, λ; r = 2, EPS0)
                    M[snlj, end - EzNo] += bz34snl(wp2, wc, kx, b1_kz2)
                end
            elseif jj == l + 1
                b2 = b2snl(p, n, l, kz, λ; EPS0)
                b4 = b4snl(p, n, l, kz, λ; EPS0)

                M[snlj, end - ExNo] += bx32snl(wp2, wc, n, kx, ξ, b2)
                M[snlj, end - EyNo] += -by32snl(wp2, wc, kx, ξ, b4)
                M[snlj, end - EzNo] += bz32snl(wp2, wc, kx, ξ, b2)
            end
        end
    end

    for idx_n in 1:(2 * N + 1), l in 1:lmax(κz)
        snl1 = _pbk_snl_index(params, s, idx_n, l, 1, N)
        M[len_sub, firstIndex + snl1] += 1
    end

    if s == length(params)
        M[len_sub, end - EzNo] += by20
    end

    return nothing
end

matrix_size(alg::BOPBK, species) = 3 * _pbk_len_sub(species, alg.N) + 6

function dispersion_matrix!(M, pb::DispersionProblem, alg::BOPBK; EPS0 = 1.0e-2, c2 = real(zero(eltype(M)) + c0)^2)
    T = _realtype(eltype(M))
    params = pb.species
    kx, kz = T(pb.kx), T(pb.kz)
    N = alg.N
    len_sub = _pbk_len_sub(params, N)
    bx10_by20 = 1im * T(ε0) * sum(sp.wp2 for sp in params)

    Mx = @view M[1:len_sub, :]
    My = @view M[(len_sub + 1):(2 * len_sub), :]
    Mz = @view M[(2 * len_sub + 1):(3 * len_sub), :]

    for s in eachindex(params)
        _pbk_add_species_Mxy!(Mx, params, s, kx, kz; MatrixNo = 1, ExNo = 5, EyNo = 4, EzNo = 3, ExyNo = 5, bx10_by20, EPS0, N, len_sub)
        _pbk_add_species_Mxy!(My, params, s, kx, kz; MatrixNo = 2, ExNo = 5, EyNo = 4, EzNo = 3, ExyNo = 4, bx10_by20, EPS0, N, len_sub)
        _pbk_add_species_Mz!(Mz, params, s, kx, kz; MatrixNo = 3, ExNo = 5, EyNo = 4, EzNo = 3, by20 = bx10_by20, EPS0, N, len_sub)
    end

    idx_Jx = len_sub
    idx_Jy = 2 * len_sub
    idx_Jz = 3 * len_sub

    inv_ε0 = 1im / T(ε0)
    M[end - 5, end - 1] += c2 * kz
    M[end - 5, idx_Jx] -= inv_ε0

    M[end - 4, end - 2] -= c2 * kz
    M[end - 4, end] += c2 * kx
    M[end - 4, idx_Jy] -= inv_ε0

    M[end - 3, end - 1] -= c2 * kx
    M[end - 3, idx_Jz] -= inv_ε0

    M[end - 2, end - 4] -= kz
    M[end - 1, end - 5] += kz
    M[end - 1, end - 3] -= kx
    M[end, end - 4] += kx

    return M
end

end
