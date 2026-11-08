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

    # End-to-end: product bi-kappa — pin a handful of alm entries.
    bk = ProductBiKappa(Electron(); T = 2555.0, kappa_perp = 200.0, kappa_para = 1.0)
    w_perp, w_para = kappa_widths(bk)
    r = HermiteVDF(bk; w_perp, w_para, N_para = 8, N_perp = 8)
    @test size(r.alm) == (9, 9)
    @test all(isfinite, r.alm)
    @test r.alm[1, 1] ≈ 1.1212568413581592 rtol = 1.0e-12
    @test r.alm[3, 1] ≈ -0.9261761839734545 rtol = 1.0e-12
    @test r.alm[1, 3] ≈ -0.005606240127199738 rtol = 1.0e-10
    @test r.alm[3, 3] ≈ 0.004630844509415841 rtol = 1.0e-10
    @test r.alm[5, 5] ≈ 0.001927003106303287 rtol = 1.0e-10
end


@testset "Hermite coefficient thresholding (atol)" begin
    # High-order expansion of a gridded Maxwellian: raw grid-quadrature noise
    # (|alm| ~ 1e-21 at l > 30) is amplified to ~1e-3 eigenvalue error;
    # atol thresholding restores the J-pole accuracy floor.
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
    e_th = hermite_expansion(fv, vz, vx, α, α; N_para = 40, N_perp = 0, atol = 1.0e-10).alm
    @test count(!iszero, e_th) == 1   # pure Maxwellian → only a₀₀ survives
    raw = hermite_expansion(fv, vz, vx, α, α; N_para = 40, N_perp = 0).alm
    @test count(!iszero, raw) > 1     # raw carries quadrature noise
    plasma = Plasma(
        Population(Electron(), HermiteVDF(; w_perp = α, w_para = α, alm = e_th); n = n_e),
        Population(Proton(), BiMaxwellian; T = T_eV, n = n_e);
        B0,
    )
    kz = 0.34 * wpe / c0_
    ωs = solve(plasma, Wavevector(0.0, kz), BOHH(; N = 2, J = 24))
    ωref = (0.09041982747948 - 0.00223551957235im) * wce   # exact Z-solver Maxwellian whistler
    @test ωs[argmin(abs.(ωs .- ωref))] ≈ ωref rtol = 1.0e-5
end
