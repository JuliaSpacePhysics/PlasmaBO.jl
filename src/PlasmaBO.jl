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

export solve
export Maxwellian
export BiKappa, BiKappa2
export gen_fv2d
export HHSolverParam
export FluidSpecies
export BranchPoint, SurfaceBranchPoint, track
export hermite_expansion
export PBK_param
export BOPBK, BOHH, BOFluid

include("types.jl")
include("utils.jl")
include("constants.jl"); using .Constants
include("formulary.jl")
include("distributions/distributions.jl")
include("Jpole.jl")
include("integral.jl")
include("hermite_expansion.jl")
include("solve.jl")
include("track.jl")

function plot_branches end

export plot_branches

end
