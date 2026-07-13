function test_hermite_coefficients_matrix(nmax)
    cHn = zeros(Float64, nmax + 1, nmax + 1)
    cHn[1, 1] = inv(sqrt(sqrt(π)))
    nmax == 0 && return cHn
    cHn[2, 2] = 2 * cHn[1, 1]
    for n in 1:(nmax - 1), k in 0:(n + 1)
        a = 2 / sqrt(n + 1)
        b = sqrt(n / (n + 1))
        cHn[n + 2, k + 1] = (k >= 1 ? a * cHn[n + 1, k] : 0.0) -
                            (k <= n - 1 ? b * cHn[n, k + 1] : 0.0)
    end
    return cHn
end

function test_hermite_a0_to_a(a0lm)
    lmax, mmax = size(a0lm)
    cHn = test_hermite_coefficients_matrix(max(lmax, mmax) - 1)
    return cHn[1:lmax, 1:lmax]' * a0lm * cHn[1:mmax, 1:mmax]
end

@testset "Hermite expansion" begin
    using PlasmaBO: hermite_H, hermite_basis
    using QuadGK: quadgk

    # H_n(x) via three-term recurrence
    @test hermite_H(0, 1.5) == 1.0
    @test hermite_H(1, 1.5) ≈ 3.0
    @test hermite_H(2, 1.5) ≈ 7.0
    @test hermite_H(3, 1.5) ≈ 9.0     # 2·1.5·7 − 2·2·3
    @test hermite_H(5, 0.7) ≈ 34.49824   # 2x·H4 − 8·H3 with x=0.7

    # ψ_l(ξ) = H_l(ξ) e^{-ξ²/2} / sqrt(2^l l! √π) is orthonormal on ℝ
    for l in 0:4, m in 0:4
        I, _ = quadgk(ξ -> hermite_basis(ξ, l) * hermite_basis(ξ, m), -12.0, 12.0; rtol = 1.0e-12)
        @test I ≈ (l == m ? 1.0 : 0.0) atol = 1.0e-10
    end

    # End-to-end: BiKappa2 — pin MATLAB power-basis values without using them in core.
    bk = BiKappa2(:e, 1.0e6, 1.0, 200.0, 2555.0; sigma = 0.0)
    g = gen_fv2d(bk)
    a0lm = hermite_expansion(g.fv, g.vz[:, 1], g.vx[1, :], g.vtz, g.vtx; Nz = 8, Nx = 8)
    alm = test_hermite_a0_to_a(a0lm)
    @test size(a0lm) == (9, 9)
    @test all(isfinite, a0lm)
    @test alm[1, 1] ≈ 1.1212568413581592 rtol = 1.0e-12
    @test alm[3, 1] ≈ -0.9261761839734545 rtol = 1.0e-12
    @test alm[1, 3] ≈ -0.005606240127199738 rtol = 1.0e-10
    @test alm[3, 3] ≈ 0.004630844509415841 rtol = 1.0e-10
    @test alm[5, 5] ≈ 0.001927003106303287 rtol = 1.0e-10
end


@testset "Hermite coefficient thresholding (atol)" begin
    using PlasmaBO: q, me
    # High-order expansion of a gridded Maxwellian: raw grid-quadrature noise
    # remains visible, but BOHH now consumes Hermite coefficients directly.
    B0 = 45.0e-9
    wce = q * B0 / me
    wpe = 100 * wce
    n_e = wpe^2 * me * 8.8541878128e-12 / q^2
    c0_ = 2.99792458e8
    α = sqrt(1.5) / 100 * c0_
    T_eV = α^2 * me / (2 * q)
    vz = collect(range(-8α, 8α, 481))
    vx = collect(range(0.0, 8α, 241))
    fv = [exp(-(z^2 + x^2) / α^2) / (π^1.5 * α^3) for z in vz, x in vx]
    e_th = hermite_expansion(fv, vz, vx, α, α; Nz = 40, Nx = 0, atol = 1.0e-10)
    @test count(!iszero, e_th) == 1   # pure Maxwellian → only a₀₀ survives
    raw = hermite_expansion(fv, vz, vx, α, α; Nz = 40, Nx = 0)
    @test count(!iszero, raw) > 1     # raw carries quadrature noise
    eM = HHSolverParam(-q, me, n_e, B0, α, α, 0.0, 0.0, e_th)
    prot = Maxwellian(:p, n_e, T_eV)
    kz = 0.34 * wpe / c0_
    ωs = solve((eM, prot), B0, 0.0, kz; N = 2, J = 24)
    ωref = (0.09041982747948 - 0.00223551957235im) * wce   # exact Z-solver Maxwellian whistler
    @test ωs[argmin(abs.(ωs .- ωref))] ≈ ωref rtol = 1.0e-5
end
