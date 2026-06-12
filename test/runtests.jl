using PlasmaBO
using Test

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(PlasmaBO)
end

include("test_track.jl")
include("test_hermite.jl")

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

@testset "Signed-density mixture (electron-deficit hole)" begin
    using PlasmaBO: q, me
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
    core = Maxwellian(:e, n_c, Tc_eV; vdz = uc / c0_)
    hole = Maxwellian(:e, -n_h, Th_eV; vdz = -vcut / c0_)   # negative density
    prot = Maxwellian(:p, n_e, Tc_eV)
    kz = 0.34 * wpe / c0_
    ωs = solve((core, hole, prot), B0, 0.0, kz; N = 2, J = 24)
    ωref = (0.0931854 + 0.00799049im) * wce
    ωbo = ωs[argmin(abs.(ωs .- ωref))]
    @test ωbo ≈ ωref rtol = 1.0e-5
    # equivalent via the direct constructor (SI velocities, negative n)
    hole_p = HHSolverParam(-q, me, -n_h, B0, 0.5α, 0.5α, -vcut, 0.0, ones(1, 1))
    core_p = HHSolverParam(-q, me, n_c, B0, α, α, uc, 0.0, ones(1, 1))
    ωs2 = solve((core_p, hole_p, prot), B0, 0.0, kz; N = 2, J = 24)
    @test ωs2[argmin(abs.(ωs2 .- ωref))] ≈ ωbo rtol = 1.0e-10
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

include("test_polarization.jl")
