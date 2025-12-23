using PlasmaBO
using Test
using Aqua

@testset "PlasmaBO.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaBO)
    end
    # Write your tests here.
end

@testset "Umeda 2012 ring beam configuration" begin
    using PlasmaBO: q, kb, ε0, me, c0

    B0 = 96.24e-9  # [Tesla]
    me_mp = 1 / 1836 # [proton mass]
    T = 51 # [eV]
    # Ring beam electrons (10% density)
    ring_beam = Species(-1.0, me_mp, 1.0e5, T; vdz = 0.1, vdr = 0.05)
    # Background electrons (90% density)
    background = Species(-1.0, me_mp, 9.0e5, T)

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
    @test filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wce, ωs) ./ wce ≈ [0.622724595851375 + 0.15724974983302017im]
end
