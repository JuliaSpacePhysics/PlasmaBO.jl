@testset "Branch tracking (track)" begin
    using PlasmaBO: DispersionSolution, BranchPoint, SurfaceBranchPoint, track

    # 1D: construct a tiny synthetic scan where the target branch is unambiguous
    ks = collect(1.0:5.0)
    ω_branch_true = ComplexF64.(10 .* ks .+ 0.1im .* ks)
    ωs_1d = map(eachindex(ks)) do i
        ωa = ω_branch_true[i]
        # Add two distractor modes far away
        [ωa, ωa + 100 + 0im, -ωa + 20im]
    end
    sol1d = DispersionSolution(ks, ωs_1d)

    k0 = ks[3]
    ω0 = ω_branch_true[3]
    k_tr, ω_tr = track(sol1d, BranchPoint(k0, ω0))
    @test ω_tr == ω_branch_true

    # 2D: build a small (k, θ) grid. Each (k, θ) has multiple candidates; the tracked
    # surface should pick the "true" mode at all points.
    θs = deg2rad.([10.0, 30.0, 50.0])
    ω_true(i, j) = ComplexF64(10 * ks[i] + 2 * j, 0.1 * ks[i] + 0.05 * j)

    ωs_2d = Matrix{Vector{ComplexF64}}(undef, length(ks), length(θs))
    for i in eachindex(ks), j in eachindex(θs)
        ωa = ω_true(i, j)
        ωs_2d[i, j] = [ωa, ωa + 80, -ωa + 10im]
    end
    sol2d = DispersionSolution(ks, θs, ωs_2d)

    seed = SurfaceBranchPoint(ks[3], θs[2], ω_true(3, 2))
    k_s, θ_s, ω_surface = track(sol2d, seed)

    ω_surface_true = [ω_true(i, j) for i in eachindex(ks), j in eachindex(θs)]
    @test ω_surface ≈ ω_surface_true
end
