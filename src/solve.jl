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
`J` controls the truncation order of the Hermite expansion used in the
solver.
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

function solve(pb::DispersionProblem, alg::BOPBK)
    M = build_pbk_dispersion_matrix(pb.species, pb.B0, pb.kx, pb.kz; N = alg.N)
    return eigvals!(M)
end

function solve(pb::DispersionProblem, alg::BOHH)
    params = HHSolverParam.(pb.species, pb.B0)
    N, J = alg.N, alg.J
    return solve_dispersion_matrix(params, pb.kx, pb.kz; N, J)
end

function solve(pb::DispersionProblem, ::BOFluid)
    M = build_fluid_dispersion_matrix(pb.species, pb.kx, pb.kz, pb.B0)
    return eigvals!(M)
end

function _ensemble_solve(f, pb)
    ks, θs = pb.ks, pb.θs
    ωs = Matrix{Vector{ComplexF64}}(undef, length(ks), length(θs))
    f(ωs)
    return DispersionSolution(ks, θs, ωs)
end

function solve(pb::EnsembleProblem, alg::BOHH)
    return _ensemble_solve(pb) do ωs
        N, J = alg.N, alg.J
        M = _zeros(pb.species, N, J)
        params = HHSolverParam.(pb.species, pb.B0)
        with_progress(pb) do ik, iθ, kx, kz
            ωs[ik, iθ] = solve_dispersion_matrix!(M, params, kx, kz; N, J)
        end
    end
end

function solve(pb::EnsembleProblem, ::BOFluid)
    return _ensemble_solve(pb) do ωs
        NN = _fluid_size(pb.species)
        M = zeros(ComplexF64, NN, NN)
        with_progress(pb) do ik, iθ, kx, kz
            fill!(M, zero(eltype(M)))
            build_fluid_dispersion_matrix!(M, pb.species, kx, kz, pb.B0)
            ωs[ik, iθ] = eigvals!(M)
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

function solve(species, B0, kx::Number, kz::Number, alg = BOHH; kw...)
    return solve(DispersionProblem(species, B0, kx, kz), alg(; kw...))
end

function solve(species, B0, ks, θs, alg = BOHH; kw...)
    return solve(EnsembleProblem(species, B0, ks, θs), alg(; kw...))
end
