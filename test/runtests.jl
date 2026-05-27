using PlasmaBO
using Test

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(PlasmaBO)
end

include("test_track.jl")

@testset "Umeda 2012 ring beam configuration" begin
    using PlasmaBO: q, me

    B0 = 96.24e-9  # [Tesla]
    T = 51 # [eV]
    ring_beam = Maxwellian(:e, 1.0e5, T; vdz = 0.1, vdr = 0.05)
    background = Maxwellian(:e, 9.0e5, T)

    species = [ring_beam, background]

    # Compute normalization
    wce = abs(B0 * q / me)
    lambdaD = Debye_length(species)
    kn = 1 / lambdaD

    # Wave vector: k*λD = 0.03, θ = 40°
    k = 0.03 / lambdaD
    θ = deg2rad(40)
    kx = k * sin(θ)
    kz = k * cos(θ)

    ωs = solve(species, B0, kx, kz; N = 6, J = 12)
    @test filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wce, ωs) ./ wce ≈ [0.6229290799953453 + 0.15687749193741884im]

    @testset "Dispersion Curve Scan" begin
        k_ranges = (0.01:0.05:0.3) .* kn
        sol = solve(species, B0, k_ranges, θ; N = 6)

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
    using PlasmaBO: qe, me, mp

    # Cold plasma parameters
    B0 = 1.0e-8  # [Tesla]
    n = 1.0e6    # [m^-3]
    T = 0.01     # [eV] - very cold

    # Wave vector
    k = 1.0  # [m^-1]
    θ = deg2rad(45)
    kx, kz = k .* sincos(θ)

    # Kinetic solver with cold Maxwellian
    e_vdf = Maxwellian(:e, n, T)
    kinetic_species = [e_vdf, Maxwellian(:p, n, T)]
    ωs_kinetic = solve(kinetic_species, B0, kx, kz; N = 2, J = 8)

    # Fluid solver
    fluid_species = [e_vdf, FluidSpecies(:p, n, T)]
    ωs_fluid = solve(fluid_species, B0, kx, kz, BOFluid)

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

@testset "Astfalk 2017 firehose instability" begin
    using PlasmaBO
    using PlasmaBO: q, me, mp
    using MAT: matopen

    B0 = 0.1
    θ = deg2rad(45)
    n = 5.0e19

    κz = 5.5
    κx = 5.5
    proton = BiKappa2(5.0e19, κz, κx, 1986.734, 993.367)
    data = gen_fv2d(proton)
    alm = hermite_expansion(data.fv, data.vz, data.vx, data.vtz, data.vtx).alm

    me_mp = 1 / 1836 # [proton mass]
    electron = Maxwellian(:e, n, 496.683)
    fpath = pkgdir(PlasmaBO, "test/firehose_Astfalk17_fvceff1.mat")
    proton_param = matopen(fpath) do file
        fvc = read(file, "fvc")
        HHSolverParam(q, mp, n, B0, fvc["vtz"], fvc["vtp"], 0.0, 0.0, fvc["alm"])
    end
    @test alm ≈ proton_param.aslm rtol = 1.0e-4

    kn = 31.0613
    k = kn / 4
    wci = proton_param.wc
    species = (proton_param, electron)
    ωs = solve(species, B0, k .* sincos(θ)...; N = 2, J = 24)
    ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
    @test imag.(ω_unstable ./ wci) ≈ [0.062373877285804444] rtol = 1.0e-3
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

    # End-to-end: BiKappa2 — pin a handful of alm entries.
    bk = BiKappa2(:e, 1.0e6, 1.0, 200.0, 2555.0; sigma = 0.0)
    g = gen_fv2d(bk)
    r = hermite_expansion(g.fv, g.vz[:, 1], g.vx[1, :], g.vtz, g.vtx; Nz = 8, Nx = 8)
    @test size(r.alm) == (9, 9)
    @test all(isfinite, r.alm)
    @test r.alm[1, 1] ≈ 1.1212568413581592 rtol = 1.0e-12
    @test r.alm[3, 1] ≈ -0.9261761839734545 rtol = 1.0e-12
    @test r.alm[1, 3] ≈ -0.005606240127199738 rtol = 1.0e-10
    @test r.alm[3, 3] ≈ 0.004630844509415841 rtol = 1.0e-10
    @test r.alm[5, 5] ≈ 0.001927003106303287 rtol = 1.0e-10
end

@testset "Arbitrary-precision solve" begin
    using DoubleFloats
    using GenericLinearAlgebra
    using LinearAlgebra
    using SpecialFunctions

    N, J = 1, 4
    sp(::Type{T}) where {T} = (HHSolverParam{T}(
        one(T), one(T), one(T), one(T), one(T), zero(T), zero(T), ones(T, 1, 1),
    ),)

    # BOHH: assembly eltype follows promoted (B0, kx, kz)
    asm(::Type{T}) where {T} = dispersion_matrix(sp(T), one(T), T(0.2), T(0.1), BOHH; N, J)
    M64 = asm(Float64)
    Md = asm(Double64)
    Mbig = setprecision(128) do; asm(BigFloat) end
    @test (eltype(M64), eltype(Md), eltype(Mbig)) ===
          (ComplexF64, Complex{Double64}, Complex{BigFloat})

    # Solve: LAPACK for F32/F64, GenericLinearAlgebra for Double64/BigFloat.
    slv(::Type{T}) where {T} = solve(sp(T), one(T), T(0.2), T(0.1), BOHH; N, J)
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
    @test eltype(solve((FluidSpecies(1.0, 1.0),), big"1", big"0.2", big"0.1", BOFluid)) === Complex{BigFloat}
    @test eltype(solve((FluidSpecies(1.0f0, 1.0f0),), 1.0f0, 0.2f0, 0.1f0, BOFluid)) === ComplexF32

    # BOPBK assembly follows promoted precision (Double64/BigFloat eigvals
    # convergence is a separate concern of the BiKappa kernel).
    bk = BiKappa2(:e, 1.0e6, 1.0, 200.0, 2555.0; sigma = 0.0)
    @test eltype(dispersion_matrix(bk, 1.0e-6, 1.0e-7, 2.0e-7, BOPBK; N = 1)) === ComplexF64
    @test eltype(dispersion_matrix(bk, Double64(1.0e-6), Double64(1.0e-7), Double64(2.0e-7), BOPBK; N = 1)) === Complex{Double64}
    @test eltype(setprecision(128) do
        dispersion_matrix(bk, big"1e-6", big"1e-7", big"2e-7", BOPBK; N = 1)
    end) === Complex{BigFloat}
end

@testset "PBK Solver" begin
    using PlasmaBO: q, me

    B0 = 1.0e-6          # [Tesla]
    θ = deg2rad(30)

    n = 2.43e6           # [m^-3]
    T = 2555.0           # [eV]
    κz = 1.0
    κx = 200.0
    electron = BiKappa2(:e, n, κz, κx, T; sigma = 0.0)

    wce = abs(B0 * q / me)
    ρce = electron.vtp / wce
    kn = 1 / ρce  # so that k/kn = kρce
    kρ_scan = range(1.0e-4, 0.3; length = 80)
    ks = kρ_scan .* kn
    sol_pbk = solve(electron, B0, ks, θ, BOPBK)
    @test length(sol_pbk.ωs) == 80
    @test maximum(real.(sol_pbk.ωs[end])) / wce ≈ 3.0667420557507734
end

@testset "Polarization and Handedness Utilities" begin
    using LinearAlgebra
    # 1. Test utilities with synthetic vectors
    # A vector where the last elements represent [Ex, Ey, Ez, Bx, By, Bz]
    # For RCP/R-wave: Ey = i*Ex => P = i (imag(P) > 0)
    v_rcp = [zeros(10)..., 1.0, 1.0im, 0.0, 0.0, 0.0, 0.0]
    @test electric_field(v_rcp) == (1.0, 1.0im, 0.0)
    @test magnetic_field(v_rcp) == (0.0, 0.0, 0.0)
    @test polarization_ratio(v_rcp) == 1.0im
    @test handedness(v_rcp) == :R

    # For LCP/L-wave: Ey = -i*Ex => P = -i (imag(P) < 0)
    v_lcp = [zeros(10)..., 1.0, -1.0im, 0.0, 0.0, 0.0, 0.0]
    @test handedness(v_lcp) == :L

    # For linear polarization: Ey = 0 => P = 0
    v_lin = [zeros(10)..., 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test handedness(v_lin) == :linear

    # 1b. Test oblique propagation (e.g. θ = 45° => kx = 1.0, kz = 1.0)
    # E_perp = Ex * cos(θ) - Ez * sin(θ) = (Ex - Ez)/√2
    # If we choose Ez = 0.0 and Ex = √2, then E_perp = 1.0.
    # If we set Ey = 1.0im, then P_transverse = Ey / E_perp = 1.0im (RCP wave).
    v_oblique = [zeros(10)..., sqrt(2), 1.0im, 0.0, 0.0, 0.0, 0.0]
    @test polarization_ratio(v_oblique, 1.0, 1.0) ≈ 1.0im
    @test polarization_ratio(v_oblique; kx = 1.0, kz = 1.0) ≈ 1.0im
    @test handedness(v_oblique, 1.0, 1.0, 1.0) == :R
    @test handedness(v_oblique, 1.0; kx = 1.0, kz = 1.0) == :R
    @test handedness(v_oblique, -1.0, 1.0, 1.0) == :L  # negative frequency reverses handedness representation
    @test handedness(v_oblique, -1.0; kx = 1.0, kz = 1.0) == :L

    # 2. Test solver with dispersion_matrix and eigen!
    using PlasmaBO: qe, me, mp
    B0 = 1.0e-4
    n = 1.0e11
    T = 1.0
    
    # Parallel propagation (along +z)
    k = 1.0
    kx, kz = 0.0, k

    e_vdf = Maxwellian(:e, n, T)
    fluid_species = [e_vdf, FluidSpecies(:p, n, T)]
    
    # Fluid solver test
    M_fluid = dispersion_matrix(fluid_species, B0, kx, kz, BOFluid)
    decomp_fluid = eigen!(M_fluid)
    @test decomp_fluid isa Eigen
    
    # Kinetic solver test
    kinetic_species = [e_vdf, Maxwellian(:p, n, T)]
    M_kinetic = dispersion_matrix(kinetic_species, B0, kx, kz, BOHH; N = 2, J = 8)
    decomp_kinetic = eigen!(M_kinetic)
    @test decomp_kinetic isa Eigen
    
    # Check that we can extract polarization from the actual solved modes
    ωs = decomp_fluid.values
    V = decomp_fluid.vectors
    
    # Let's find any mode with a non-zero real frequency
    idx = findfirst(ω -> isfinite(ω) && abs(real(ω)) > 1.0e3, ωs)
    if idx !== nothing
        v = V[:, idx]
        Ex, Ey, Ez = electric_field(v)
        @test isfinite(polarization_ratio(v))
        @test handedness(v) in (:R, :L, :linear)
    end

    # Helper to construct a minimal eigenvector with specified E and B field components
    function make_vector(Ex, Ey, Ez; Bx=0.0, By=0.0, Bz=0.0)
        v = zeros(ComplexF64, 9)
        v[4] = Ex; v[5] = Ey; v[6] = Ez
        v[7] = Bx; v[8] = By; v[9] = Bz
        return v
    end

    v_rcp = make_vector(1.0 + 0im, 0.5im, 0.0)
    v_lcp = make_vector(1.0 + 0im, -0.5im, 0.0)

    # polarization_ratio: k==0 fallback returns Ey/Ex
    @test polarization_ratio(make_vector(2.0 + 0im, 6.0 + 0im, 0.0), 0.0, 0.0) ≈ 3.0

    # handedness with k==0: E_perp = Ex branch (covers L102)
    @test handedness(v_rcp, 1.0, 0.0, 0.0) == :R

    # handedness: negligible E fields → :linear (covers L112)
    @test handedness(make_vector(0.0 + 0im, 0.0 + 0im, 0.0), 1.0, 0.0, 1.0) == :linear

    # handedness: R and L polarizations
    @test handedness(v_rcp, 1.0, 0.0, 1.0) == :R
    @test handedness(v_lcp, 1.0, 0.0, 1.0) == :L

    # handedness: significant real-valued P (|P| above threshold) with zero ellipticity → :linear (covers L129)
    @test handedness(make_vector(1.0 + 0im, 1.0 + 0im, 0.0), 1.0, 0.0, 1.0) == :linear

    # Dispatch: handedness(v; kx, kz) forwards correctly (covers L135/L137)
    @test handedness(v_rcp) == handedness(v_rcp, 1.0, 0.0, 1.0)
    @test handedness(v_rcp; kx=0.0, kz=1.0) == handedness(v_rcp, 1.0, 0.0, 1.0)

    # Dispatch: handedness(v, ω) defaults to parallel propagation (covers L143)
    ω = 2.0 + 0im
    @test handedness(v_rcp, ω) == handedness(v_rcp, ω, 0.0, 1.0)
    # Dispatch: handedness(v, ω; kx, kz) with explicit kwargs (covers L145)
    @test handedness(v_rcp, ω; kx=0.0, kz=1.0) == handedness(v_rcp, ω, 0.0, 1.0)
end
