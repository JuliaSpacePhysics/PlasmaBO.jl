# PlasmaBO

[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/refs/heads/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("PlasmaBO")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://JuliaSpacePhysics.github.io/PlasmaBO.jl/dev/)


## Features and Roadmap

- Hermite-Hermite (HH) expansion for arbitrary/analytic distributions
    - [x] Maxwellian / BiMaxwellian
    - [ ] Kappa / BiKappa / product Bikappa
- [ ] Generalized plasma dispersion function (GPDF)
- [ ] Integration with [VelocityDistributionFunctions.jl](https://github.com/JuliaSpacePhysics/VelocityDistributionFunctions.jl) and observation / simulation data
- [ ] Multi-fluid solver
- [ ] Reformulate as a `LinearProblem` for use with `SciML` [LinearSolve](https://docs.sciml.ai/LinearSolve/stable/)
- [ ] Faster eigenvalue solver using Krylov methods ([Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) / [KrylovKit](https://github.com/Jutho/KrylovKit.jl))

## Elsewhere

- [hsxie/BO-Arbitrary](https://github.com/hsxie/boarbitrary/tree/main): Extension of the kinetic electromagnetic magnetized dispersion relation solver [PDRK](https://github.com/hsxie/pdrk)/[BO](https://github.com/hsxie/bo) to arbitrary distributions.
- [danielver02/ALPS](https://github.com/danielver02/ALPS): The Arbitrary Linear Plasma Solver that solves the Vlasov-Maxwell dispersion relation in hot (even relativistic) magnetised plasma.
- [pastfalk/LEOPARD](https://github.com/pastfalk/LEOPARD): Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions
