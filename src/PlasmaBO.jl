module PlasmaBO

using Gamma: gamma, loggamma
using QuadGK: quadgk, quadgk!
using Bumper: @no_escape, @alloc
using LinearAlgebra
using Tullio: @tullio
using ProgressMeter: @showprogress, next!
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
export electric_field, magnetic_field, polarization_ratio, handedness

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
