```@meta
CurrentModule = PlasmaBO
```

# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)
[![version](https://juliahub.com/docs/General/PlasmaBO/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlasmaBO)

Documentation for [PlasmaBO](https://github.com/JuliaSpacePhysics/PlasmaBO.jl).

## Installation

```julia
using Pkg
Pkg.add("PlasmaBO")
```

## Features

- Hermite-Hermite (BO-HH) expansion solver for arbitrary/analytic distributions
    - Maxwellian / BiMaxwellian
- Analytic PBK (BO-PBK) solver for kappa distributions
- Multi-fluid solver
- Integration with [`ChargedParticles.jl`](https://juliaplasma.github.io/ChargedParticles.jl/dev/)

## Usage Examples

The matrix eigenvalue method ([xieEfficientFrameworkSolving2025](@citet), [xiePDRKGeneralKinetic2016](@citet)) finds all wave modes simultaneously by transforming the dispersion relation into a matrix eigenvalue problem using J-pole approximation for the plasma dispersion function.

This approach is more efficient to find multiple modes at once, and doesn't require initial guesses for the root finder.

Check out the [ring beam instability example](ringbeam_Umeda12.md) for detailed usage instructions, also see [firehose instability example](firehose_Astfalk17.md) for using with arbitrary velocity distributions, [BO-PBK example](rlp_Cattaert07.md) for using with kappa distributions (BO-PBK), [cold plasma example](cold_plasma.md) for comparing kinetic and fluid solvers, and [dispersion surface tracking example](dispersion_surface.md) for 2D scanning and mode tracking.

### Solvers

BO-PBK ([`BOPBK`](@ref)) is an analytic, distribution-aware eigen-solver optimized for kappa plasmas, whereas BO-Arbitrary ([`BOHH`](@ref)) is a universal but numerically heavier framework that approximates any distribution at the cost of efficiency and low-κ accuracy.

```@docs; canonical = false
BOPBK
BOHH
BOFluid
```

## Notation & Assumptions

The formulation (code) is valid for non-relativistic, arbitrary gyrotropic distributions.

- **Coordinates**
  - **`z`**: direction parallel to the background magnetic field (**`B0`**).
  - **`x`**: one perpendicular direction (any perpendicular direction is equivalent).
- **Velocities**
  - **`vz`, `vx`**: particle velocity components along `z` and `x`.
  - **`vdz`, `vdx`**: drift/bulk velocity components along `z` and `x` (when present in a distribution parameterization).
  - **`vtz`, `vtx`**: thermal speeds along `z` and `x`.
- **Wave vector**
  - **`θ`**: propagation angle between `k` and `B0`.
  - **`k∥ = k cos(θ)`**, **`k⊥ = k sin(θ)`**.

## References

[xieRapidComputationPlasma2024](@citet), [xieBO20Plasma2021](@citet), [xiePDRFGeneralDispersion2014](@citet), [xieGeneralizedPlasmaDispersion2013](@citet), 

```@bibliography
```

## API Reference

```@index
```

```@autodocs
Modules = [PlasmaBO]
```
