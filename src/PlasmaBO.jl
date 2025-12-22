module PlasmaBO

using SpecialFunctions
using FFTW
using QuadGK
using LinearAlgebra

export Species, solve_kinetic_dispersion
export FluidSpecies, FluidSolverParams, solve_fluid_dispersion
export BranchPoint, track_dispersion_branch, track_dispersion_branches

include("types.jl")
include("utils.jl")
include("constants.jl")
include("formulary.jl")
include("Jpole.jl")
include("integral.jl")
include("matrix_solver.jl")
include("fluid_solver.jl")
include("eigenvalue_filtering.jl")

end
