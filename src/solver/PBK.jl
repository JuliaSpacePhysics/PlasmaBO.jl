module PBK
using QuadGK: quadgk
using SpecialFunctions: loggamma, gamma, besselj
using ..Constants

function PBK_param(s, B0)
    T = Float64
    u0 = T(s.vdz) * c0
    sigma = T(s.sigma)
    κz = Int(s.κz)
    κx = T(s.κx)
    m = s.m
    q = s.q
    wp = abs(q) * sqrt(s.n / (m * ε0))
    wc = (T(B0) * q) / m

    vtz = sqrt((2 - 1 / κz) * qe * s.Tz / m)
    vtx = sqrt((2 - 2 / κx) * qe * s.Tp / m / (1 + sigma))
    rhoc = vtx / abs(wc)

    return (; u0, sigma, κz, κx, wp, wc, vtz, vtx, rhoc)
end

function _pbk_species_list(species)
    species isa AbstractVector && return species
    species isa Tuple && return species
    return (species,)
end

_pbk_per_n(κ) = (κ + 1) * (κ + 4) ÷ 2

function _pbk_len_total(species, Ns)
    species = _pbk_species_list(species)
    len = 0
    for sp in species
        len += (2 * Ns + 1) * _pbk_per_n(Int(sp.κz))
    end
    return len
end

_pbk_len_sub(species, N) = _pbk_len_total(species, N) + 1

function _pbk_species_offset(species, s, Ns)
    species = _pbk_species_list(species)
    off = 0
    for i in 1:(s - 1)
        sp = species[i]
        off += (2 * Ns + 1) * _pbk_per_n(Int(sp.κz))
    end
    return off
end

function _pbk_snl_index(species, s, idx_n, l, j, Ns)
    species = _pbk_species_list(species)
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

function _csn(p, n, kz)
    κz = p.κz
    return n * p.wc + kz * p.u0 - 1im * sqrt(κz) * kz * p.vtz
end

function _log_csl(κz, l, kz, vtz)
    return loggamma(κz + 1) + loggamma(2 * κz - l + 2) -
        loggamma(κz - l + 2) - loggamma(2 * κz + 1) +
        (l - 1) * log(2im * sqrt(κz) * kz * vtz)
end

function bsnl(p, n, l, kz)
    log_csl = _log_csl(p.κz, l, kz, p.vtz)
    return -n * p.wc * exp(log_csl)
end

function bsl(p, l, kz)
    log_csl = _log_csl(p.κz, l, kz, p.vtz)
    return -0.5 * l * kz^2 * p.vtx^2 * exp(log_csl)
end

function funcS_pbk(p, n, λ, p1, p2, p3, p4; EPS0 = 1.0e-2)
    κx = p.κx
    sgm = p.sigma
    eps0 = eps(Float64)

    if λ > EPS0
        S0 = 4.0 * (2 * λ)^(-sgm - 2) * κx^(-sgm - 1) * exp(loggamma(κx + sgm + 1) - loggamma(κx)) / gamma(sgm + 1)
        function integrand1(x)
            jn = besselj(n, x)
            djn = 0.5 * (besselj(n - 1, x) - besselj(n + 1, x))
            num = (x + eps0)^(2 * sgm + p3) * jn^p1 * djn^p2
            den = (1 + 0.5 * x^2 / λ / κx)^(κx + sgm + p4)
            return num / den
        end

        tmp = quadgk(integrand1, 0.0, Inf)[1]
        return S0 * tmp
    else
        S0 = 4.0 * (2 * λ)^(p3 / 2 - 1.5) * κx^(-sgm - 1) * exp(loggamma(κx + sgm + 1) - loggamma(κx)) / gamma(sgm + 1)
        a = sqrt(2 * λ)
        function integrand2(x)
            jn = besselj(n, x * a)
            djn = 0.5 * (besselj(n - 1, x * a) - besselj(n + 1, x * a))
            num = (x + eps0)^(2 * sgm + p3) * jn^p1 * djn^p2
            den = (1 + x^2 / κx)^(κx + sgm + p4)
            return num / den
        end
        tmp = quadgk(integrand2, 0.0, Inf)[1]
        return S0 * tmp
    end
end

function b1snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S1 = funcS_pbk(p, n, λ, 2, 0, 1, 2; EPS0)
    S7 = funcS_pbk(p, n, λ, 2, 0, -1, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S17 = (κx + sgm + 1) / κx * S1 - 2 * sgm * λ * S7
    return S17 * bsnl(p, n, l, kz)
end

function b2snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S2 = funcS_pbk(p, n, λ, 2, 0, 1, 1; EPS0)
    return S2 * bsl(p, l, kz)
end

function b3snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S3 = funcS_pbk(p, n, λ, 1, 1, 2, 2; EPS0)
    S8 = funcS_pbk(p, n, λ, 1, 1, 0, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S38 = (κx + sgm + 1) / κx * S3 - 2 * sgm * λ * S8
    return S38 * bsnl(p, n, l, kz)
end

function b4snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S4 = funcS_pbk(p, n, λ, 1, 1, 2, 1; EPS0)
    return S4 * bsl(p, l, kz)
end

function b5snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S5 = funcS_pbk(p, n, λ, 0, 2, 3, 2; EPS0)
    S9 = funcS_pbk(p, n, λ, 0, 2, 1, 1; EPS0)
    κx = p.κx
    sgm = p.sigma
    S59 = (κx + sgm + 1) / κx * S5 - 2 * sgm * λ * S9
    return S59 * bsnl(p, n, l, kz)
end

function b6snl(p, n, l, kz, λ; EPS0 = 1.0e-2)
    S6 = funcS_pbk(p, n, λ, 0, 2, 3, 1; EPS0)
    return S6 * bsl(p, l, kz)
end

bx11snl(wp, n, b1) = -1im * ε0 * wp^2 * n^2 * b1
bx12snl(wp, n, b2) = -1im * ε0 * wp^2 * n^2 * b2
bx21snl(wp, n, b3) = ε0 * wp^2 * n * b3
bx22snl(wp, n, b4) = ε0 * wp^2 * n * b4
bx31snl(wp, wc, n, tanθ, csn, b1, b2) =
    -1im * ε0 * tanθ * wp^2 * n * (b1 * (csn - n * wc) + b2) / wc
bx32snl(wp, wc, n, tanθ, csn, b2) =
    -1im * ε0 * tanθ * wp^2 * n * b2 * (csn - n * wc) / wc
bx33snl(wp, wc, n, tanθ, b1) =
    -1im * ε0 * tanθ * wp^2 * n * b1 / wc

by21snl(wp, b5) = -1im * ε0 * wp^2 * b5
by22snl(wp, b6) = -1im * ε0 * wp^2 * b6
by31snl(wp, wc, n, tanθ, csn, b3, b4) =
    -ε0 * tanθ * wp^2 * (b3 * (csn - n * wc) + b4) / wc
by32snl(wp, wc, n, tanθ, csn, b4) =
    -ε0 * tanθ * wp^2 * b4 * (csn - n * wc) / wc
by33snl(wp, wc, n, tanθ, b3) =
    -ε0 * tanθ * wp^2 * b3 / wc

bz31snl(wp, wc, n, tanθ, csn, b1, b2) =
    -1im * ε0 * tanθ^2 * wp^2 * (b1 * (csn - n * wc)^2 + 2 * b2 * (csn - n * wc)) / wc^2
bz32snl(wp, wc, n, tanθ, csn, b2) =
    -1im * ε0 * tanθ^2 * wp^2 * b2 * (csn - n * wc)^2 / wc^2
bz33snl(wp, wc, n, tanθ, csn, b1, b2) =
    -1im * ε0 * tanθ^2 * wp^2 * (2 * b1 * (csn - n * wc) + b2) / wc^2
bz34snl(wp, wc, tanθ, b1) =
    -1im * ε0 * tanθ^2 * wp^2 * b1 / wc^2

lmax(κz; threshold = 10) = min(κz + 1, threshold)

function _pbk_add_species_Mxy!(
        M, params, s,
        kx, kz;
        MatrixNo, ExNo, EyNo, EzNo, ExyNo,
        bx10_by20, EPS0, len_sub, N
    )
    p = params[s]
    firstIndex = _pbk_block_first_index(len_sub, MatrixNo)
    λ = 0.5 * kx^2 * (p.rhoc^2)
    κz = p.κz

    wp = p.wp
    wc = p.wc
    tanθ = kx / kz

    for (idx_n, n) in enumerate(-N:N)
        csn = _csn(p, n, kz)
        for l in 1:lmax(κz)
            for jj in 1:(l + 1)
                snlj = _pbk_snl_index(params, s, idx_n, l, jj, N)

                M[snlj, firstIndex + snlj] += csn
                if jj < l + 1
                    M[snlj, firstIndex + snlj + 1] += 1
                end

                if jj == l
                    if MatrixNo == 1
                        b1 = b1snl(p, n, l, kz, λ; EPS0)
                        b2 = b2snl(p, n, l, kz, λ; EPS0)
                        b3 = b3snl(p, n, l, kz, λ; EPS0)

                        M[snlj, end - ExNo] += bx11snl(wp, n, b1)
                        M[snlj, end - EyNo] += bx21snl(wp, n, b3)
                        M[snlj, end - EzNo] += bx31snl(wp, wc, n, tanθ, csn, b1, b2)

                        if l <= κz
                            b1_lp1 = b1snl(p, n, l + 1, kz, λ; EPS0)
                            M[snlj, end - EzNo] += bx33snl(wp, wc, n, tanθ, b1_lp1)
                        end
                    else
                        b3 = b3snl(p, n, l, kz, λ; EPS0)
                        b4 = b4snl(p, n, l, kz, λ; EPS0)
                        b5 = b5snl(p, n, l, kz, λ; EPS0)

                        M[snlj, end - ExNo] += -bx21snl(wp, n, b3)
                        M[snlj, end - EyNo] += by21snl(wp, b5)
                        M[snlj, end - EzNo] += by31snl(wp, wc, n, tanθ, csn, b3, b4)

                        if l <= κz
                            b3_lp1 = b3snl(p, n, l + 1, kz, λ; EPS0)
                            M[snlj, end - EzNo] += by33snl(wp, wc, n, tanθ, b3_lp1)
                        end
                    end
                elseif jj == l + 1
                    if MatrixNo == 1
                        b2 = b2snl(p, n, l, kz, λ; EPS0)
                        b4 = b4snl(p, n, l, kz, λ; EPS0)

                        M[snlj, end - ExNo] += bx12snl(wp, n, b2)
                        M[snlj, end - EyNo] += bx22snl(wp, n, b4)
                        M[snlj, end - EzNo] += bx32snl(wp, wc, n, tanθ, csn, b2)
                    else
                        b4 = b4snl(p, n, l, kz, λ; EPS0)
                        b6 = b6snl(p, n, l, kz, λ; EPS0)

                        M[snlj, end - ExNo] += -bx22snl(wp, n, b4)
                        M[snlj, end - EyNo] += by22snl(wp, b6)
                        M[snlj, end - EzNo] += by32snl(wp, wc, n, tanθ, csn, b4)
                    end
                end
            end
        end
    end

    for idx_n in 1:(2 * N + 1)
        for l in 1:(κz + 1)
            snl1 = _pbk_snl_index(params, s, idx_n, l, 1, N)
            M[len_sub, firstIndex + snl1] += 1
        end
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
    λ = 0.5 * kx^2 * (p.rhoc^2)
    κz = p.κz
    wp = p.wp
    wc = p.wc
    tanθ = kx / kz
    for (idx_n, n) in enumerate(-N:N)
        csn = _csn(p, n, kz)
        for l in 1:lmax(κz)
            for jj in 1:(l + 1)
                snlj = _pbk_snl_index(params, s, idx_n, l, jj, N)

                M[snlj, firstIndex + snlj] += csn
                if jj < l + 1
                    M[snlj, firstIndex + snlj + 1] += 1
                end

                if jj == l
                    b1 = b1snl(p, n, l, kz, λ; EPS0)
                    b2 = b2snl(p, n, l, kz, λ; EPS0)
                    b3 = b3snl(p, n, l, kz, λ; EPS0)
                    b4 = b4snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += bx31snl(wp, wc, n, tanθ, csn, b1, b2)
                    M[snlj, end - EyNo] += -by31snl(wp, wc, n, tanθ, csn, b3, b4)
                    M[snlj, end - EzNo] += bz31snl(wp, wc, n, tanθ, csn, b1, b2)

                    if l <= κz
                        b1_lp1 = b1snl(p, n, l + 1, kz, λ; EPS0)
                        b2_lp1 = b2snl(p, n, l + 1, kz, λ; EPS0)
                        b3_lp1 = b3snl(p, n, l + 1, kz, λ; EPS0)

                        M[snlj, end - ExNo] += bx33snl(wp, wc, n, tanθ, b1_lp1)
                        M[snlj, end - EyNo] += -by33snl(wp, wc, n, tanθ, b3_lp1)
                        M[snlj, end - EzNo] += bz33snl(wp, wc, n, tanθ, csn, b1_lp1, b2_lp1)
                    end

                    if l <= κz - 1
                        b1_lp2 = b1snl(p, n, l + 2, kz, λ; EPS0)
                        M[snlj, end - EzNo] += bz34snl(wp, wc, tanθ, b1_lp2)
                    end
                elseif jj == l + 1
                    b2 = b2snl(p, n, l, kz, λ; EPS0)
                    b4 = b4snl(p, n, l, kz, λ; EPS0)

                    M[snlj, end - ExNo] += bx32snl(wp, wc, n, tanθ, csn, b2)
                    M[snlj, end - EyNo] += -by32snl(wp, wc, n, tanθ, csn, b4)
                    M[snlj, end - EzNo] += bz32snl(wp, wc, n, tanθ, csn, b2)
                end
            end
        end
    end

    for idx_n in 1:(2 * N + 1)
        for l in 1:lmax(κz)
            snl1 = _pbk_snl_index(params, s, idx_n, l, 1, N)
            M[len_sub, firstIndex + snl1] += 1
        end
    end

    if s == length(params)
        M[len_sub, end - EzNo] += by20
    end

    return nothing
end

function _pbk_Mxy(params, kx, kz; N, kw...)
    len_sub = _pbk_len_sub(params, N)
    len_col = 3 * len_sub + 6
    M = zeros(ComplexF64, len_sub, len_col)
    for s in eachindex(params)
        _pbk_add_species_Mxy!(
            M,
            params, s, kx, kz;
            len_sub, N, kw...
        )
    end
    return M
end

function _pbk_Mz(params, kx, kz; N, kw...)
    len_sub = _pbk_len_sub(params, N)
    len_col = 3 * len_sub + 6
    M = zeros(ComplexF64, len_sub, len_col)

    for s in eachindex(params)
        _pbk_add_species_Mz!(M, params, s, kx, kz; len_sub, N, kw...)
    end

    return M
end

function build_pbk_dispersion_matrix(species, B0, kx, kz; N = 2, EPS0 = 1.0e-2, c2 = c0^2)
    species = _pbk_species_list(species)
    len_sub = _pbk_len_sub(species, N)
    params = PBK_param.(species, B0)
    wp2_sum = sum(sp.wp^2 for sp in params)

    bx10_by20 = 1im * ε0 * wp2_sum
    by20 = bx10_by20

    Mx = _pbk_Mxy(params, kx, kz; MatrixNo = 1, ExNo = 5, EyNo = 4, EzNo = 3, ExyNo = 5, bx10_by20, EPS0, N)
    My = _pbk_Mxy(params, kx, kz; MatrixNo = 2, ExNo = 5, EyNo = 4, EzNo = 3, ExyNo = 4, bx10_by20, EPS0, N)
    Mz = _pbk_Mz(params, kx, kz; MatrixNo = 3, ExNo = 5, EyNo = 4, EzNo = 3, by20, EPS0, N)

    O = zeros(ComplexF64, 6, size(Mz, 2))
    M = vcat(Mx, My, Mz, O)

    idx_Jx = len_sub
    idx_Jy = 2 * len_sub
    idx_Jz = 3 * len_sub

    M[end - 5, end - 1] += c2 * kz
    M[end - 5, idx_Jx] -= 1im / ε0

    M[end - 4, end - 2] -= c2 * kz
    M[end - 4, end] += c2 * kx
    M[end - 4, idx_Jy] -= 1im / ε0

    M[end - 3, end - 1] -= c2 * kx
    M[end - 3, idx_Jz] -= 1im / ε0

    M[end - 2, end - 4] -= kz
    M[end - 1, end - 5] += kz
    M[end - 1, end - 3] -= kx
    M[end, end - 4] += kx

    return M
end
end
