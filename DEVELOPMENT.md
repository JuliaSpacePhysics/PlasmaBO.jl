[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl)

## Roadmap

- [ ] Extend BO-Product-Bi-Kappa (PBK) solver to non-integer κ (For integer κ, $Z_κ(ξ)$ has an exact finite closed-form expansion).
- [ ] Generalized plasma dispersion function (GPDF)
- [ ] Integration with [VelocityDistributionFunctions.jl](https://github.com/JuliaSpacePhysics/VelocityDistributionFunctions.jl) and observation / simulation data  
- [ ] Faster eigenvalue solver using Krylov methods ([Arpack](https://github.com/JuliaLinearAlgebra/Arpack.jl) / [KrylovKit](https://github.com/Jutho/KrylovKit.jl), ref: [Eigen solvers](https://docs.sciml.ai/BifurcationKit/stable/eigensolver/))
- [ ] GPU Acceleration / Parallelization / Sparse matrix optimizations
- [ ] Reformulate as a `SciMLProblem` for use with `SciML` (ref: [LinearSolve](https://docs.sciml.ai/LinearSolve/stable/), [ApproxFun.jl](https://juliaapproximation.github.io/ApproxFun.jl/stable/generated/Eigenvalue/))
- [ ] Relativistic support
    - [ ] Derive the dielectric/susceptibility tensor from the relativistic Vlasov–Maxwell equation expressed in momentum space;
    - [ ] relativistic analogs of plasma dispersion function;
    - [ ] Choose basis functions in momentum space (relativistic Maxwell–Jüttner weighted bases).
- [ ] Better handling of long-tailed distributions