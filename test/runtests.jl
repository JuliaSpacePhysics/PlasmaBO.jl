using PlasmaBO
using PlasmaBO: hermite_expansion
using PlasmaPopulations: gyrofrequency, thermal_speed, debye_length, vth_perp_para, kappa_widths, pdf
using PlasmaPopulations: E_SI as q, ME_SI as me, MP_SI as mp, C_SI as c0, EPS0_SI as ε0
using Test

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(PlasmaBO)
end

include("test_track.jl")
include("test_hermite.jl")

@testset "Umeda 2012 ring beam configuration" begin

    B0 = 96.24e-9  # [Tesla]
    T = 51 # [eV]
    e = Electron()
    vth = thermal_speed(e, T)
    plasma = Plasma(
        Population(e, GaussianRing(; vr_perp = 0.05 * c0, w_perp = vth, vth_para = vth, vd_para = 0.1 * c0); n = 1.0e5),
        Population(e, BiMaxwellian(; vth); n = 9.0e5);
        B0,
    )

    # Compute normalization
    wce = abs(B0 * q / me)
    lambdaD = debye_length(e, 1.0e6, T)
    kn = 1 / lambdaD

    # Wave vector: k*λD = 0.03, θ = 40°
    k = 0.03 / lambdaD
    θ = deg2rad(40)
    kx = k * sin(θ)
    kz = k * cos(θ)

    ωs = solve(plasma, Wavevector(kx, kz), BOHH(; N = 6, J = 12))
    @test filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wce, ωs) ./ wce ≈ [0.6229290799953453 + 0.15687749193741884im]

    @testset "Dispersion Curve Scan" begin
        k_ranges = (0.01:0.05:0.3) .* kn
        sol = solve(plasma, AngleSweep(; k = k_ranges, theta = θ), BOHH(; N = 6))

        initial_points = [
            (0.1 * kn, 0.3im * wce),
            (0.1 * kn, 0.1im * wce),
            (0.2 * kn, 0.25im * wce),
        ]

        branches = track.(sol, initial_points)
        ω_maxs = [maximum(imag.(b[2])) for b in branches]
        @test ω_maxs ≈ [5425.78665660744, 1740.4141552857307, 5425.78665660744]
    end
end

@testset "Fluid vs Kinetic: cold plasma limit" begin

    # Cold plasma parameters
    B0 = 1.0e-8  # [Tesla]
    n = 1.0e6    # [m^-3]
    T = 0.01     # [eV] - very cold

    # Wave vector
    k = 1.0  # [m^-1]
    θ = deg2rad(45)
    kx, kz = k .* sincos(θ)

    plasma = Plasma(Population(Electron(), BiMaxwellian; T, n), Population(Proton(), BiMaxwellian; T, n); B0)
    ωs_kinetic = solve(plasma, Wavevector(kx, kz), BOHH(; N = 2, J = 8))
    ωs_fluid = solve(plasma, Wavevector(kx, kz), BOFluid())

    # Compare electromagnetic wave modes (highest frequency modes, near ±c*k)
    # These should agree well in cold plasma limit
    filter_em = ω -> isfinite(ω) && abs(real(ω)) > 1.0e8
    fluid_em = sort(filter(filter_em, ωs_fluid), by = ω -> real(ω))
    kinetic_em = sort(filter(filter_em, ωs_kinetic), by = ω -> real(ω))

    # Should have same number of EM modes (light waves)
    @test length(fluid_em) == length(kinetic_em) == 4

    # EM wave frequencies should match within 1e-6 relative error
    @test fluid_em[1:2] ≈ kinetic_em[1:2] rtol = 1.0e-6
end

@testset "BOFluid: warm drifting species, parallel longitudinal modes" begin
    B0, n, k = 1.0e-6, 1.0e12, 1.0e-1
    e, p = Electron(), Proton()
    plasma = Plasma(
        Population(e, BiMaxwellian(e; T = 100.0, vd_para = 1.0e6); n),
        Population(p, BiMaxwellian(p; T = 10.0, vd_para = -1.0e5); n); B0
    )
    γ = (3.0, 1.0)
    ωs = solve(plasma, Wavevector(0.0, k), BOFluid(; gamma_para = γ))
    # 1 = Σ ωp² / ((ω - k vd)² - γ k² T/m), T = m vth²/2
    terms(ω) = map(plasma.populations, γ) do s, γs
        (; q, m) = s.particle
        s.n * q^2 / (ε0 * m) / ((ω - k * s.vdf.vd_para)^2 - γs * k^2 * vth_perp_para(s.vdf)[2]^2 / 2)
    end
    residual(ω) = abs(1 - sum(terms(ω))) / (1 + sum(abs, terms(ω)))
    @test count(ω -> residual(ω) < 1.0e-8, ωs) == 4
    @test_throws ArgumentError solve(plasma, Wavevector(0.0, k), BOFluid(; gamma_para = [3.0]))
end

@testset "Astfalk 2017 firehose instability" begin
    using PlasmaBO
    using MAT: matopen

    B0 = 0.1
    θ = deg2rad(45)
    n = 5.0e19

    p = Proton()
    proton = ProductBiKappa(p; T_perp = 993.367, T_para = 1986.734, kappa_perp = 5.5, kappa_para = 5.5)
    w_perp, w_para = kappa_widths(proton)
    hermite = HermiteVDF(proton; w_perp, w_para)

    fpath = pkgdir(PlasmaBO, "test/firehose_Astfalk17_fvceff1.mat")
    ref = matopen(fpath) do file
        fvc = read(file, "fvc")
        HermiteVDF(; w_perp = fvc["vtp"], w_para = fvc["vtz"], alm = fvc["alm"])
    end
    @test hermite.alm ≈ ref.alm rtol = 1.0e-4

    kn = 31.0613
    k = kn / 4
    wci = gyrofrequency(p, B0)
    plasma = Plasma(Population(p, ref; n), Population(Electron(), BiMaxwellian; T = 496.683, n); B0)
    ωs = solve(plasma, Wavevector(k .* sincos(θ)...), BOHH(; N = 2, J = 24))
    ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
    @test imag.(ω_unstable ./ wci) ≈ [0.062373877285804444] rtol = 1.0e-3
end

@testset "Signed-density mixture (electron-deficit hole)" begin
    # Core Maxwellian minus a drifting-Maxwellian hole at v∥ = −v_cut: the
    # sunward-deficit whistler instability.
    # Reference values from an exact plasma-dispersion-function (Z) solver
    # for the same signed components.
    # Setup: ω_pe/|Ω_e| = 100, β_c = 1.5, κ = v_cut/α_c = 2, hole width 0.5α_c,
    # full depth.
    B0 = 45.0e-9
    wce = q * B0 / me
    wpe = 100 * wce
    n_e = wpe^2 * me * 8.8541878128e-12 / q^2
    c0_ = 2.99792458e8
    η = 0.5^3 * exp(-4.0)
    n_c = n_e / (1 - η)
    n_h = η * n_c
    α = sqrt(1.5 * n_e / n_c) / 100 * c0_
    vcut = 2.0 * α
    uc = -η * vcut
    Tc_eV = α^2 * me / (2 * q)
    Th_eV = (0.5 * α)^2 * me / (2 * q)
    e, p = Electron(), Proton()
    plasma = Plasma(
        Population(e, BiMaxwellian; T = Tc_eV, vd_para = uc, n = n_c),
        Population(e, BiMaxwellian; T = Th_eV, vd_para = -vcut, n = -n_h),   # negative density
        Population(p, BiMaxwellian; T = Tc_eV, n = n_e);
        B0,
    )
    kz = 0.34 * wpe / c0_
    ωs = solve(plasma, Wavevector(0.0, kz), BOHH(; N = 2, J = 24))
    ωref = (0.0931854 + 0.00799049im) * wce
    @test ωs[argmin(abs.(ωs .- ωref))] ≈ ωref rtol = 1.0e-5
end

@testset "B0 given as a vector" begin
    sp = (Population(Proton(), BiMaxwellian(; vth = 1.0e5); n = 5.0e6), Population(Electron(), BiMaxwellian(; vth = 4.0e6); n = 5.0e6))
    M(B0) = dispersion_matrix(Plasma(sp...; B0), Wavevector(1.0e-6, 2.0e-6), BOHH())
    @test M([3.0e-9, 0.0, 4.0e-9]) ≈ M(5.0e-9)
end

@testset "Arbitrary-precision solve" begin
    using DoubleFloats
    using GenericLinearAlgebra
    using LinearAlgebra
    using SpecialFunctions

    N, J = 1, 4
    # unit gyro- and plasma frequency, unit thermal speeds
    unit(::Type{T}) where {T} = Plasma(
        Population(Particle(one(T), one(T)), HermiteVDF(; w_perp = one(T), w_para = one(T), alm = ones(T, 1, 1)); n = T(ε0));
        B0 = one(T),
    )

    # BOHH: assembly eltype follows promoted (B0, kx, kz)
    asm(::Type{T}) where {T} = dispersion_matrix(unit(T), Wavevector(T(0.2), T(0.1)), BOHH(; N, J))
    M64 = asm(Float64)
    Md = asm(Double64)
    Mbig = setprecision(128) do; asm(BigFloat) end
    @test (eltype(M64), eltype(Md), eltype(Mbig)) ===
          (ComplexF64, Complex{Double64}, Complex{BigFloat})

    # Solve: LAPACK for F32/F64, GenericLinearAlgebra for Double64/BigFloat.
    slv(::Type{T}) where {T} = solve(unit(T), Wavevector(T(0.2), T(0.1)), BOHH(; N, J))
    ω32 = slv(Float32)
    ωd = slv(Double64)
    ωbig = setprecision(128) do; slv(BigFloat) end
    @test eltype(ω32) === ComplexF32 && all(isfinite, ω32)
    @test eltype(ωd) === Complex{Double64} && all(isfinite, ωd)
    @test eltype(ωbig) === Complex{BigFloat} && all(isfinite, ωbig)

    # Susceptibility block (the 3x3 species sum on the field columns) converges
    # to the BigFloat reference much faster in extended precision.
    SNJ = (2N + 1) * J
    SNJ1 = SNJ + 1
    block(M) = M[[SNJ + 1, SNJ + SNJ1 + 1, SNJ + 2SNJ1 + 1], (3SNJ1 + 1):(3SNJ1 + 3)]
    Bref = setprecision(256) do; block(asm(BigFloat)) end
    err(M) = norm(Complex{BigFloat}.(block(M)) - Bref) / norm(Bref)
    @test err(Md) < err(M64) * big"1e-10"

    # BOFluid: same generic path; precision follows wave-parameter promotion.
    fluid(::Type{T}, B0) where {T} = Plasma(
        Population(Particle(T(q), T(mp)), BiMaxwellian(; vth = T(1.0e4)); n = one(T));
        B0,
    )
    @test eltype(solve(fluid(Float64, big"1"), Wavevector(big"0.2", big"0.1"), BOFluid())) === Complex{BigFloat}
    @test eltype(solve(fluid(Float32, 1.0f0), Wavevector(0.2f0, 0.1f0), BOFluid())) === ComplexF32

    # BOPBK assembly follows promoted precision (Double64/BigFloat eigvals
    # convergence is a separate concern of the BiKappa kernel).
    e = Electron()
    bk(B0) = Plasma(Population(e, ProductBiKappa; T = 2555.0, kappa_perp = 200.0, kappa_para = 1.0, n = 1.0e6); B0)
    @test eltype(dispersion_matrix(bk(1.0e-6), Wavevector(1.0e-7, 2.0e-7), BOPBK(; N = 1))) === ComplexF64
    @test eltype(dispersion_matrix(bk(Double64(1.0e-6)), Wavevector(Double64(1.0e-7), Double64(2.0e-7)), BOPBK(; N = 1))) === Complex{Double64}
    @test eltype(setprecision(128) do
        dispersion_matrix(bk(big"1e-6"), Wavevector(big"1e-7", big"2e-7"), BOPBK(; N = 1))
    end) === Complex{BigFloat}
end

@testset "PBK Solver" begin

    B0 = 1.0e-6          # [Tesla]
    θ = deg2rad(30)

    n = 2.43e6           # [m^-3]
    T = 2555.0           # [eV]
    e = Electron()
    electron = ProductBiKappa(e; T, kappa_perp = 200.0, kappa_para = 1.0)

    wce = abs(B0 * q / me)
    ρce = kappa_widths(electron)[1] / wce
    kn = 1 / ρce  # so that k/kn = kρce
    kρ_scan = range(1.0e-4, 0.3; length = 80)
    ks = kρ_scan .* kn
    sol_pbk = solve(Plasma(Population(e, electron; n); B0), AngleSweep(; k = ks, theta = θ), BOPBK())
    @test length(sol_pbk.omega) == 80
    @test maximum(real.(sol_pbk.omega[end])) / wce ≈ 3.0667420557507734
end

@testset "Negative density subtracts a population" begin
    e, B0, n = Electron(), 1.0e-6, 2.43e6
    T = 2555.0
    kx, kz = 1.0e-3, 2.0e-3
    for (vdf, alg) in (
            (ProductBiKappa(e; T, kappa_perp = 200.0, kappa_para = 1.0), BOPBK),
            (BiMaxwellian(e; T), BOFluid),
            (BiMaxwellian(e; T), BOHH),
        )
        ref = solve(Plasma(Population(e, vdf; n); B0), Wavevector(kx, kz), alg())
        split = solve(Plasma(Population(e, vdf; n = 2n), Population(e, vdf; n = -n); B0), Wavevector(kx, kz), alg())
        @test all(ω -> minimum(abs, split .- ω) ≤ 1.0e-5 * maximum(abs, ref), ref)
    end
end

@testset "A wavevector's form is explicit" begin
    e = Electron()
    plasma = Plasma(Population(e, BiMaxwellian; T = 100.0, n = 1.0e6); B0 = 1.0e-6)
    k, θ = 1.0e-3, deg2rad(30)
    alg = BOHH()
    cartesian = solve(plasma, Wavevector(k .* sincos(θ)...), alg)
    polar = solve(plasma, AngleSweep(; k = [k], theta = θ), alg).omega[1, 1]
    @test all(ω -> minimum(abs, polar .- ω) ≤ 1.0e-6 * maximum(abs, cartesian), cartesian)
    @test_throws MethodError solve(plasma, k, θ, alg)
    @test_throws MethodError solve(plasma, Wavevector(k, k))
end

@testset "BOFluid takes a ring's pressure from its second moment" begin
    ring = GaussianRing(; vr_perp = 1.0e6, w_perp = 1.0e6, vth_para = 1.0e6)
    vth_perp, vth_para = vth_perp_para(ring)
    M(d) = dispersion_matrix(Plasma(Population(Electron(), d; n = 1.0e6); B0 = 1.0e-6), Wavevector(1.0e-3, 2.0e-3), BOFluid())
    @test M(ring) ≈ M(BiMaxwellian(; vth_perp, vth_para))
end

@testset "Exactly perpendicular propagation" begin
    e, p = Electron(), Proton()
    B0, n = 1.0e-5, 1.0e12
    warm = Plasma(Population(e, BiMaxwellian(e; T = 1.0, vd_para = 1.0e4); n), Population(p, BiMaxwellian(p; T = 1.0); n); B0)
    kappa = Plasma(Population(e, ProductBiKappa(e; T = 2555.0, kappa_perp = 4.5, kappa_para = 3, loss_cone = 1.0, vd_para = 1.0e6); n = 2.43e6); B0 = 1.0e-6)
    kx = 0.1

    # The matrix is continuous into k∥ = 0
    for (plasma, alg) in ((warm, BOHH()), (kappa, BOPBK(; N = 1)))
        M0 = dispersion_matrix(plasma, Wavevector(kx, 0.0), alg)
        M1 = dispersion_matrix(plasma, Wavevector(kx, 1.0e-10 * kx), alg)
        @test all(isfinite, M0)
        @test all(isapprox.(M0, M1; rtol = 1.0e-7, atol = 1.0e-10 * maximum(abs, M0)))
    end

    # O-mode: ω² = ωp² + c²k² for kρ ≪ 1
    ωp2 = n * q^2 / ε0 * (1 / me + 1 / mp)
    ωs = solve(warm, Wavevector(kx, 0.0), BOHH())
    @test minimum(abs.(ωs .- sqrt(ωp2 + (c0 * kx)^2))) / sqrt(ωp2) < 1.0e-3

    # Undamped modes at kρ = 1.5, ωpe/ωce ≈ 2; the references are VlasovMaxwellDispersion.jl's roots [ω/ωce]
    B0, n, T = 1.0e-6, 3.9e7, 100.0
    wce = abs(gyrofrequency(e, B0))
    k = Wavevector(1.5 * wce / thermal_speed(e, T), 0.0)
    for (vdf, alg, roots) in (
            (BiMaxwellian, BOHH(; N = 8), (1.34938, 2.40623, 3.18098)),
            (ProductBiKappa, BOPBK(; N = 8), (1.39918, 2.39876, 3.15771)),
        )
        kw = vdf === BiMaxwellian ? (;) : (; kappa_perp = 4.0, kappa_para = 4)
        plasma = Plasma(Population(e, vdf; T, n, kw...), Population(p, vdf; T, n, kw...); B0)
        ωs = solve(plasma, k, alg) ./ wce
        @test all(ω -> minimum(abs.(ωs .- ω)) < 3.0e-4, roots)
    end
end

@testset "GridVDF through BOHH" begin
    B0, n = 0.1, 5.0e19
    p, e = Proton(), Electron()
    k = Wavevector((31.0613 / 4) .* sincos(deg2rad(45))...)
    wci = gyrofrequency(p, B0)
    table(d, (w_perp, w_para), vd = 0.0) = let v_perp = range(0, 12w_perp, 241), v_para = range(vd - 12w_para, vd + 12w_para, 481)
        GridVDF(v_perp, v_para, [pdf(d, x, z) for x in v_perp, z in v_para])
    end
    growth(proton, electron) =
        maximum(imag, filter(isfinite, solve(Plasma(Population(p, proton; n), Population(e, electron; n); B0), k, BOHH(; N = 2, J = 24)))) / wci

    # A tabulated drifting bi-Maxwellian against the same model solved exactly
    electron = BiMaxwellian(e; T = 496.683)
    maxwellian = BiMaxwellian(p; T_perp = 993.367, T_para = 1986.734, vd_para = 1.0e5)
    @test growth(table(maxwellian, vth_perp_para(maxwellian), 1.0e5), electron) ≈ growth(maxwellian, electron) rtol = 1.0e-2

    # A tabulated kappa expands like its model; the reference is VlasovMaxwellDispersion.jl's susceptibility
    kappa = ProductBiKappa(p; T_perp = 993.367, T_para = 1986.734, kappa_perp = 6.0, kappa_para = 6)
    tabulated = growth(table(kappa, vth_perp_para(kappa)), electron)
    @test tabulated ≈ growth(kappa, electron) rtol = 1.0e-2
    @test tabulated ≈ 0.06314 rtol = 0.02
end

include("test_polarization.jl")
