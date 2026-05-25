module PlasmaBO

using Gamma: gamma, loggamma
using QuadGK: quadgk, quadgk!
using Bumper: @no_escape, @alloc
using LinearAlgebra
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
export BOPBK, BOHH, BOFluid
export electric_field, magnetic_field, polarization_ratio, handedness, dispersion_matrix

function dispersion_matrix end
function dispersion_matrix! end
function matrix_size end
prepare(alg, species, B0) = species

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
include("relativistic.jl")
using .RelativisticLongitudinal

export RelativisticLongitudinal
export solve_relativistic_direct, solve_relativistic_matrix, solve_maxwellian_matrix
export juttner_dielectric, langmuir_root, omega0_asymptotic

function plot_branches end

export plot_branches

end
