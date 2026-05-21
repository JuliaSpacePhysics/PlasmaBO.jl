include("solver/PBK.jl"); using .PBK: build_pbk_dispersion_matrix, PBK_param
include("solver/HH.jl")
include("solver/fluid.jl")

"""
    BOPBK(; N = 2)

Dispersion solver using the PBK matrix formulation.

`N` controls the truncation order of the cyclotron harmonic index used to build
the dispersion matrix.
"""
@kwdef struct BOPBK <: AbstractDispersionAlgorithm
    N::Int = 2
end

"""
    BOHH(; N = 2, J = 8)

Dispersion solver using the Hermite-Hankel (HH) matrix formulation.

`N` controls the truncation order of the cyclotron harmonic index.
`J` controls the truncation order of the number of poles for Z-function approximation

This method transforms the dispersion relation into a matrix eigenvalue problem
using J-pole approximation for the plasma dispersion function, allowing
simultaneous computation of all wave modes.
"""
@kwdef struct BOHH <: AbstractDispersionAlgorithm
    N::Int = 2
    J::Int = 8
end

"""
    BOFluid()

Dispersion solver using the multi-fluid electromagnetic matrix formulation.
"""
@kwdef struct BOFluid <: AbstractDispersionAlgorithm end

function dispersion_matrix(pb::DispersionProblem, alg::BOPBK)
    return build_pbk_dispersion_matrix(pb.species, pb.B0, pb.kx, pb.kz; N = alg.N)
end

function solve(pb::DispersionProblem, alg::BOPBK)
    M = dispersion_matrix(pb, alg)
    return eigvals!(M)
end

zero!(M) = fill!(M, zero(eltype(M)))

function dispersion_matrix(pb::DispersionProblem, alg::BOHH)
    params = HHSolverParam.(pb.species, pb.B0)
    N, J = alg.N, alg.J
    n = _size(length(pb.species), N, J)
    M = zeros(ComplexF64, n, n)
    build_dispersion_matrix!(M, params, pb.kx, pb.kz; N, J)
    return M
end

function solve(pb::DispersionProblem, alg::BOHH)
    params = HHSolverParam.(pb.species, pb.B0)
    N, J = alg.N, alg.J
    n = _size(length(pb.species), N, J)
    return @no_escape begin
        M = @alloc(ComplexF64, n, n)
        build_dispersion_matrix!(zero!(M), params, pb.kx, pb.kz; N, J)
        eigvals!(M)
    end
end

function dispersion_matrix(pb::DispersionProblem, ::BOFluid)
    return build_fluid_dispersion_matrix(pb.species, pb.kx, pb.kz, pb.B0)
end

function solve(pb::DispersionProblem, alg::BOFluid)
    M = dispersion_matrix(pb, alg)
    return eigvals!(M)
end

function _ensemble_solve(f, pb)
    ks, θs = pb.ks, pb.θs
    ωs = Matrix{Vector{ComplexF64}}(undef, length(ks), length(θs))
    solve_with_threads(1) do
        f(ωs)
    end
    return DispersionSolution(ks, θs, ωs)
end

function solve(pb::EnsembleProblem, alg::BOHH)
    return _ensemble_solve(pb) do ωs
        N, J = alg.N, alg.J
        params = HHSolverParam.(pb.species, pb.B0)
        n = _size(length(pb.species), N, J)
        with_progress(pb) do ik, iθ, kx, kz
            @no_escape begin
                M = @alloc(ComplexF64, n, n)
                build_dispersion_matrix!(zero!(M), params, kx, kz; N, J)
                ωs[ik, iθ] = eigvals!(M)
            end
        end
    end
end

function solve(pb::EnsembleProblem, ::BOFluid)
    return _ensemble_solve(pb) do ωs
        NN = _fluid_size(pb.species)
        with_progress(pb) do ik, iθ, kx, kz
            @no_escape begin
                M = @alloc(ComplexF64, NN, NN)
                build_fluid_dispersion_matrix!(zero!(M), pb.species, kx, kz, pb.B0)
                ωs[ik, iθ] = eigvals!(M)
            end
        end
    end
end

function solve(pb::EnsembleProblem, alg::BOPBK)
    return _ensemble_solve(pb) do ωs
        with_progress(pb) do ik, iθ, kx, kz
            spb = DispersionProblem(pb.species, pb.B0, kx, kz)
            ωs[ik, iθ] = solve(spb, alg)
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

function solve(species, B0, ks, θs, alg = BOHH; kw...)
    return solve(EnsembleProblem(species, B0, ks, θs), alg(; kw...))
end
