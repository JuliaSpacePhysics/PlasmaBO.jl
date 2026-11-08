module PlasmaBO

using QuadGK: quadgk, quadgk!
using Bumper: @no_escape, @alloc
import CommonSolve: solve
using LinearAlgebra
using ProgressMeter: @showprogress, next!
using PlasmaPopulations
using PlasmaPopulations: si_velocity, C_SI as c0, EPS0_SI as ε0

export Particle, Electron, Positron, Proton
export Population, Plasma
export BiMaxwellian, GaussianRing, BiKappa, ProductBiKappa, Cold, GridVDF
export Wavevector, AngleSweep
export solve, dispersion_matrix
export BOPBK, BOHH, BOFluid
export HermiteVDF
export track, plot_branches
export polarization_ratio, handedness

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public hermite_expansion, BranchPoint, SurfaceBranchPoint, DispersionSolution"))

function dispersion_matrix end
function dispersion_matrix! end
function matrix_size end
function prepare end
function plot_branches end

include("types.jl")
include("utils.jl")
include("vdfs.jl")
include("Jpole.jl")
include("special.jl")
include("integral.jl")
include("hermite_expansion.jl")
include("solve.jl")
include("track.jl")

end
