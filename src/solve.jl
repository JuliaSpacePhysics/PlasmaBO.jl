include("solver/PBK.jl"); using .PBK: BOPBK
include("solver/HH.jl")
include("solver/fluid.jl")

zero!(M) = fill!(M, zero(eltype(M)))

function dispersion_matrix(pb::DispersionProblem, alg; kw...)
    n = matrix_size(alg, pb.species)
    return dispersion_matrix!(zeros(ComplexF64, n, n), pb, alg; kw...)
end

function solve(pb::DispersionProblem, alg)
    n = matrix_size(alg, pb.species)
    return @no_escape begin
        M = @alloc(ComplexF64, n, n)
        dispersion_matrix!(zero!(M), pb, alg)
        eigvals!(M)
    end
end

function _ensemble_solve(f, pb)
    ks, θs = pb.ks, pb.θs
    ωs = Matrix{Vector{ComplexF64}}(undef, length(ks), length(θs))
    solve_with_threads(1) do
        f(ωs)
    end
    return DispersionSolution(ks, θs, ωs)
end

function solve(pb::EnsembleProblem, alg;
        progress = false,
        progress_name = "Solving dispersion (k, θ)...",
        progress_steps = 4)
    return _ensemble_solve(pb) do ωs
        species = prepare(alg, pb.species, pb.B0)
        n = matrix_size(alg, species)
        with_progress(pb; progress, name = progress_name, progress_steps) do ik, iθ, kx, kz
            @no_escape begin
                M = @alloc(ComplexF64, n, n)
                sub = DispersionProblem(species, pb.B0, kx, kz)
                dispersion_matrix!(zero!(M), sub, alg)
                ωs[ik, iθ] = eigvals!(M)
            end
        end
    end
end

"""
    solve(species, B0, kx, kz, alg = BOHH; kw...)
    solve(species, B0, ks, θs, alg = BOHH; kw...)

Convenience interface for solving a single wavevector `(kx, kz)` or
    dispersion scan over multiple wavevectors `(ks, θs)`.

Keyword arguments are passed to the algorithm constructor. Defaults to `BOHH` solver.

See also: [`BOHH`](@ref), [`BOPBK`](@ref), [`BOFluid`](@ref)
"""
function solve end

function dispersion_matrix(species, B0, kx::Number, kz::Number, alg = BOHH; kw...)
    return dispersion_matrix(DispersionProblem(species, B0, kx, kz), alg(; kw...))
end

function solve(species, B0, kx::Number, kz::Number, alg = BOHH; kw...)
    return solve(DispersionProblem(species, B0, kx, kz), alg(; kw...))
end

function solve(species, B0, ks, θs, alg = BOHH;
        progress = false,
        progress_name = "Solving dispersion (k, θ)...",
        progress_steps = 1,
        kw...)
    return solve(EnsembleProblem(species, B0, ks, θs), alg(; kw...);
        progress, progress_name, progress_steps)
end
