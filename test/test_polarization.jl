@testset "Polarization and Handedness Utilities" begin
    using LinearAlgebra
    electric_field(v) = (v[end - 5], v[end - 4], v[end - 3])
    magnetic_field(v) = (v[end - 2], v[end - 1], v[end])
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
    function make_vector(Ex, Ey, Ez; Bx = 0.0, By = 0.0, Bz = 0.0)
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
    @test handedness(v_rcp; kx = 0.0, kz = 1.0) == handedness(v_rcp, 1.0, 0.0, 1.0)

    # Dispatch: handedness(v, ω) defaults to parallel propagation (covers L143)
    ω = 2.0 + 0im
    @test handedness(v_rcp, ω) == handedness(v_rcp, ω, 0.0, 1.0)
    # Dispatch: handedness(v, ω; kx, kz) with explicit kwargs (covers L145)
    @test handedness(v_rcp, ω; kx = 0.0, kz = 1.0) == handedness(v_rcp, ω, 0.0, 1.0)
end
