module PlasmaBO

using QuadGK: quadgk, quadgk!
using Bumper: @no_escape, @alloc
using LinearAlgebra
using PlasmaPopulations
using PlasmaPopulations: si_velocity
using Reexport: @reexport
@reexport using PlasmaDispersionBase
using PlasmaDispersionBase: normalized, realtype
import PlasmaDispersionBase: pointsolver

export dispersion_matrix
export BOPBK, BOHH, BOFluid
export HermiteVDF
export polarization_ratio, handedness

function dispersion_matrix end
function dispersion_matrix! end
function matrix_size end
function prepare end

include("utils.jl")
include("vdfs.jl")
include("Jpole.jl")
include("special.jl")
include("integral.jl")
include("hermite_expansion.jl")
include("solve.jl")

end
