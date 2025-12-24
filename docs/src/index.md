```@meta
CurrentModule = PlasmaBO
```

# PlasmaBO

Documentation for [PlasmaBO](https://github.com/JuliaSpacePhysics/PlasmaBO.jl).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaSpacePhysics/PlasmaBO.jl")
```

## Features

- [x] Hermite-Hermite (HH) expansion for arbitrary/analytic distributions
    - [x] Maxwellian / BiMaxwellian
    - [x] Kappa / BiKappa / product Bikappa
- [x] Integration with [`ChargedParticles.jl`](https://juliaplasma.github.io/ChargedParticles.jl/dev/)

## Usage Examples

The matrix eigenvalue method finds all wave modes simultaneously by transforming the dispersion relation into a matrix eigenvalue problem using J-pole approximation for the plasma dispersion function.

This approach is more efficient to find multiple modes at once, and doesn't require initial guesses for the root finder.

Check out the [ring beam instability example](ringbeam_Umeda12.md) for detailed usage instructions, also see [firehose instability example](firehose_Astfalk17.md) for using with arbitrary velocity distributions.

## References

[xieEfficientFrameworkSolving2025](@citet), [xieRapidComputationPlasma2024](@citet), [xieBO20Plasma2021](@citet), [xiePDRKGeneralKinetic2016](@citet), [xiePDRFGeneralDispersion2014](@citet), [xieGeneralizedPlasmaDispersion2013](@citet), 

```@bibliography
```

## API Reference

```@index
```

```@autodocs
Modules = [PlasmaBO]
```
