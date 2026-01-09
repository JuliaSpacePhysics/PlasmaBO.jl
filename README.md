# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)
[![version](https://juliahub.com/docs/General/PlasmaBO/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlasmaBO)

[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/refs/heads/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("PlasmaBO")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://JuliaSpacePhysics.github.io/PlasmaBO.jl/dev/)


## Features and Roadmap

- [x] Hermite-Hermite (HH) basis solver for arbitrary/analytic distributions
    - [x] Maxwellian / BiMaxwellian
    - [x] Kappa / BiKappa / product Bikappa
- [x] BO-Product-Bi-Kappa (PBK) solver for kappa distributions
    - [ ] Extend to non-integer κ (For integer κ, $Z_κ(ξ)$ has an exact finite closed-form expansion).
- [ ] Generalized plasma dispersion function (GPDF)
- [ ] Integration with [VelocityDistributionFunctions.jl](https://github.com/JuliaSpacePhysics/VelocityDistributionFunctions.jl) and observation / simulation data
- [x] Multi-fluid solver
- [ ] Faster eigenvalue solver using Krylov methods ([Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) / [KrylovKit](https://github.com/Jutho/KrylovKit.jl), ref: [Eigen solvers](https://docs.sciml.ai/BifurcationKit/stable/eigensolver/))
- [ ] GPU Acceleration / Parallelization / Sparse matrix optimizations
- [ ] Reformulate as a `SciMLProblem` for use with `SciML` (ref: [LinearSolve](https://docs.sciml.ai/LinearSolve/stable/), [ApproxFun.jl](https://juliaapproximation.github.io/ApproxFun.jl/stable/generated/Eigenvalue/))
- [ ] Relativistic support
    - [ ] Rederive the dielectric/susceptibility tensor from the relativistic Vlasov–Maxwell equation expressed in momentum space;
    - [ ] relativistic analogs of plasma dispersion function;
    - [ ] Choose basis functions in momentum space (relativistic Maxwell–Jüttner weighted bases).
- [ ] Better handling of long-tailed distributions

## Elsewhere

- [hsxie/BO-Arbitrary](https://github.com/hsxie/boarbitrary/tree/main): Extension of the kinetic electromagnetic magnetized dispersion relation solver [PDRK](https://github.com/hsxie/pdrk)/[BO](https://github.com/hsxie/bo) to arbitrary distributions (MATLAB)
    - [liangwang0734/xenon](https://github.com/liangwang0734/xenon): A matriX-based dispErsioN relatiON solver (Python). Very limited in functionality (only electrostatic for kinetic BiMaxwellian plasmas).
- [LinearMaxwellVlasov.jl](https://github.com/jwscook/LinearMaxwellVlasov.jl): A Julia implementation of linear Maxwell-Vlasov solvers
- [danielver02/ALPS](https://github.com/danielver02/ALPS): The Arbitrary Linear Plasma Solver that solves the Vlasov-Maxwell dispersion relation in hot (even relativistic) magnetised plasma (Fortran)
- [pastfalk/LEOPARD](https://github.com/pastfalk/LEOPARD): Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions (Fortran)
- [Drakicy/MPDES](https://github.com/Drakicy/MPDES): Magnetized Plasma Dispersion Equation Solver (MPDES) (MATLAB)
