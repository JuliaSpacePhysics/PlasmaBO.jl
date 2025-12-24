module PlasmaBO

using SpecialFunctions
using SpecialFunctions: gamma, erf, erfc
using QuadGK: quadgk
using LinearAlgebra
using Tullio: @tullio
using ProgressMeter: @showprogress
using ChargedParticles: charge, mass, charge_number, mass_number, particle, ParticleLike
import ChargedParticles as CP
using Unitful
# using KrylovKit
using ArnoldiMethod

export solve_kinetic_dispersion
export Maxwellian
export BiKappa, BiKappa2
export gen_fv2d
export HHSolverParam
export FluidSpecies, FluidSolverParams, solve_fluid_dispersion
export BranchPoint, track_dispersion_branch, track_dispersion_branches
export hermite_expansion

include("types.jl")
include("utils.jl")
include("constants.jl")
include("formulary.jl")
include("distributions/distributions.jl")
include("Jpole.jl")
include("integral.jl")
include("hermite_expansion.jl")
include("matrix_solver.jl")
include("fluid_solver.jl")
include("eigenvalue_filtering.jl")

end
