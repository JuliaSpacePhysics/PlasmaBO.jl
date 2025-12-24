using PlasmaBO
using Test

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(PlasmaBO)
end

@testset "Umeda 2012 ring beam configuration" begin
    using PlasmaBO: q, kb, ε0, me, c0

    B0 = 96.24e-9  # [Tesla]
    me_mp = 1 / 1836 # [proton mass]
    T = 51 # [eV]
    # Ring beam electrons (10% density)
    ring_beam = Maxwellian(-1.0, me_mp, 1.0e5, T; vdz = 0.1, vdr = 0.05)
    # Background electrons (90% density)
    background = Maxwellian(-1.0, me_mp, 9.0e5, T)

    species = [ring_beam, background]

    # Compute normalization
    wce = abs(B0 * q / me)
    lambdaD = Debye_length(species)

    # Wave vector: k*λD = 0.03, θ = 40°
    k = 0.03 / lambdaD
    θ = deg2rad(40)
    kx = k * sin(θ)
    kz = k * cos(θ)

    # Solve using matrix eigenvalue method
    # J=12 provides good accuracy (J-pole approximation order)
    ωs = solve_kinetic_dispersion(species, B0, kx, kz; N = 6, J = 12)
    @test filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wce, ωs) ./ wce ≈ [0.6229225524157314 + 0.15686544364848592im]
end


@testset "Astfalk 2017 firehose instability" begin
    using PlasmaBO
    using PlasmaBO: q, kb, ε0, me, c0, mp
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
    electron = Maxwellian(-1.0, me_mp, n, 496.683)
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
    ωs = solve_kinetic_dispersion(species, B0, k .* sincos(θ)...; N = 2, J = 24)
    ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
    @test imag.(ω_unstable ./ wci) ≈ [0.062373877285804444] rtol = 1.0e-3
end
