using PlasmaBO
using Test

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(PlasmaBO)
end

@testset "Umeda 2012 ring beam configuration" begin
    using PlasmaBO: q, me

    B0 = 96.24e-9  # [Tesla]
    T = 51 # [eV]
    # Ring beam electrons (10% density)
    ring_beam = Maxwellian(:e, 1.0e5, T; vdz = 0.1, vdr = 0.05)
    # Background electrons (90% density)
    background = Maxwellian(:e, 9.0e5, T)

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
    @test filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wce, ωs) ./ wce ≈ [0.6229290799953453 + 0.15687749193741884im]
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
    ωs_kinetic = solve_kinetic_dispersion(kinetic_species, B0, kx, kz; N = 2, J = 8)

    # Fluid solver
    fluid_species = [e_vdf, FluidSpecies(:p, n, T)]
    ωs_fluid = solve_fluid_dispersion(fluid_species, B0, kx, kz)

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
    ωs = solve_kinetic_dispersion(species, B0, k .* sincos(θ)...; N = 2, J = 24)
    ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
    @test imag.(ω_unstable ./ wci) ≈ [0.062373877285804444] rtol = 1.0e-3
end
