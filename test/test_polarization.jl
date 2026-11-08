@testset "Polarization and Handedness Utilities" begin
    using LinearAlgebra
    field_vector(Ex, Ey, Ez) = ComplexF64[0, 0, 0, Ex, Ey, Ez, 0, 0, 0]
    parallel, oblique = Wavevector(0.0, 1.0), Wavevector(1.0, 1.0)

    # Ey = ±i Ex: P = ±i (imag(P) > 0 is R)
    v_R, v_L = field_vector(1, 1im, 0), field_vector(1, -1im, 0)
    @test polarization_ratio(v_R, parallel) == 1.0im
    @test handedness(v_R, 1.0, parallel) == :R
    @test handedness(v_L, 1.0, parallel) == :L
    @test handedness(field_vector(1, 0, 0), 1.0, parallel) == :linear
    @test handedness(field_vector(1, 1, 0), 1.0, parallel) == :linear   # real P, zero ellipticity
    @test handedness(field_vector(0, 0, 0), 1.0, parallel) == :linear   # no field

    # θ = 45°: E⊥ = (Ex - Ez)/√2, so Ex = √2 gives E⊥ = 1; the sign of Re ω reverses the sense.
    v_oblique = field_vector(sqrt(2), 1im, 0)
    @test polarization_ratio(v_oblique, oblique) ≈ 1.0im
    @test handedness(v_oblique, 1.0, oblique) == :R
    @test handedness(v_oblique, -1.0, oblique) == :L

    # k = 0 falls back to Ey/Ex
    @test polarization_ratio(field_vector(2, 6, 0), Wavevector(0.0, 0.0)) ≈ 3.0
    @test handedness(v_R, 1.0, Wavevector(0.0, 0.0)) == :R

    # Eigenvectors of a solver's matrix end in the fields
    plasma = Plasma(Population(Electron(), BiMaxwellian; T = 1.0, n = 1.0e11), Population(Proton(), BiMaxwellian; T = 1.0, n = 1.0e11); B0 = 1.0e-4)
    k = Wavevector(0.0, 1.0)
    for alg in (BOFluid(), BOHH(; N = 2, J = 8))
        decomp = eigen!(dispersion_matrix(plasma, k, alg))
        i = findfirst(ω -> isfinite(ω) && abs(real(ω)) > 1.0e3, decomp.values)
        @test isfinite(polarization_ratio(decomp.vectors[:, i], k))
        @test handedness(decomp.vectors[:, i], decomp.values[i], k) in (:R, :L, :linear)
    end
end
