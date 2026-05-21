# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)
[![version](https://juliahub.com/docs/General/PlasmaBO/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlasmaBO)

[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaBO.jl)

## Quickstart

`solve` returns all wave modes at once as eigenvalues of a matrix dispersion relation — no initial guess required. Define species, pick a solver (`BOHH` kinetic (Hermite-Hermite basis for arbitrary distributions), `BOPBK` product Bi-Kappa, `BOFluid` multi-fluid), and call:

```julia
using Pkg; Pkg.add("PlasmaBO")
using PlasmaBO

B0 = 100e-9                       # background field [T]
n, T = 8.7e6, 2.857e-3            # density [m^-3], temperature [eV]
species = (Maxwellian(:p, n, T), Maxwellian(:e, n, T))

# Single (kx, kz): returns Vector{ComplexF64} of all ω roots
kx, kz = 0.0, 1e-6
ωs = solve(species, B0, kx, kz)            # default: BOHH kinetic
ωs_fluid = solve(species, B0, kx, kz, BOFluid)

# Scan over |k| and propagation angle θ
ks = 10 .^ range(-7, -4, length = 50)
θs = deg2rad.(0:15:90)
sol = solve(species, B0, ks, θs, BOHH)     # sol.ωs[ik, iθ] :: Vector{ComplexF64}
```

`ω = ωr + iγ`: `γ > 0` is an instability, `γ < 0` is damping. Use `Maxwellian` / `BiKappa` for kinetic species; `ChargedParticles.jl` shorthands like `Maxwellian("O-18 3+", n, T)` work for arbitrary ions.

See [documentation](https://JuliaSpacePhysics.github.io/PlasmaBO.jl/dev/) for worked examples: cold plasma (kinetic vs fluid), ring-beam, firehose, ion-beam, and dispersion-surface tracking.

## Elsewhere

- [hsxie/BO-Arbitrary](https://github.com/hsxie/boarbitrary/tree/main): Extension of the kinetic electromagnetic magnetized dispersion relation solver [PDRK](https://github.com/hsxie/pdrk)/[BO](https://github.com/hsxie/bo) to arbitrary distributions (MATLAB)
  - [liangwang0734/xenon](https://github.com/liangwang0734/xenon): A matriX-based dispErsioN relatiON solver (Python). Very limited in functionality (only electrostatic for kinetic BiMaxwellian plasmas).
- [LinearMaxwellVlasov.jl](https://github.com/jwscook/LinearMaxwellVlasov.jl): A Julia implementation of linear Maxwell-Vlasov solvers
- [danielver02/ALPS](https://github.com/danielver02/ALPS): The Arbitrary Linear Plasma Solver that solves the Vlasov-Maxwell dispersion relation in hot (even relativistic) magnetised plasma (Fortran)
- [pastfalk/LEOPARD](https://github.com/pastfalk/LEOPARD): Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions (Fortran)
- [Drakicy/MPDES](https://github.com/Drakicy/MPDES): Magnetized Plasma Dispersion Equation Solver (MPDES) (MATLAB)
