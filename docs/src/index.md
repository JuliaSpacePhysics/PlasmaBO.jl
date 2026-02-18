```@meta
CurrentModule = PlasmaBO
```

# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)
[![version](https://juliahub.com/docs/General/PlasmaBO/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlasmaBO)

## Installation

```julia
using Pkg
Pkg.add("PlasmaBO")
```

## Problem Definition

The core problem is solving the linear dispersion relation for waves in a magnetized plasma:

```math
\det(\mathbf{D}(\omega, \mathbf{k})) = 0
```
where $\mathbf{D}$ is the dispersion tensor, $\omega$ is the complex frequency, and $\mathbf{k}$ is the wave vector.

**Input**:
- **Plasma Parameters**: For each species $s$, density $n_s$, temperature $T_s$ (anisotropic $T_{\parallel}, T_{\perp}$), drift velocity $v_{d,s}$, charge $q_s$, and mass $m_s$.
    - **Note**: The code formulation is valid for non-relativistic, arbitrary gyrotropic distributions.
- **Coordinates & Fields**:
    - **`z`**: direction parallel to the background magnetic field (**`B0`**).
    - **`x`**: one perpendicular direction (any perpendicular direction is equivalent).
    - ** Velocities**: `vz` (parallel), `vx` (perpendicular), thermal speeds `vtz`, `vtx`.
- **Wave Parameters**: Wave vector $\mathbf{k} = (k_\perp, 0, k_\parallel)$.
    - **`θ`**: propagation angle between `k` and `B0`.
    - **`k∥ = k cos(θ)`**, **`k⊥ = k sin(θ)`**.

**Output**:
- **Complex Frequency**: $\omega = \omega_r + i\gamma$.
    - $\omega_r$: Real frequency (oscillation).
    - $\gamma$: Growth rate (instability if $\gamma > 0$, damping if $\gamma < 0$).

## Methodology

Solving the kinetic dispersion relation is traditionally challenging because:
1.  **Transcendental Equations**: The dispersion tensor involves the Plasma Dispersion Function $Z(\zeta)$, making the equation $D(\omega) = 0$ transcendental and highly nonlinear.
2.  **Root Finding Difficulty**: Traditional iterative methods (like Newton-Raphson) require good initial guesses. Without them, it's easy to miss modes or converge to the wrong root.
3.  **Multimode Nature**: Plasmas support multiple wave modes simultaneously (e.g., Alfvén, Whistlers, Langmuir), making it hard to track them all.

**PlasmaBO.jl** overcomes these challenges using a **Matrix Eigenvalue Method** ([xieEfficientFrameworkSolving2025](@citet), [xiePDRKGeneralKinetic2016](@citet)):

1.  **Rational Approximation**: We approximate the Plasma Dispersion Function $Z(\zeta)$ using a multi-pole (J-pole) expansion:
    ```math
    Z(\zeta) \approx \sum_{j=1}^J \frac{b_j}{\zeta - c_j}
    ```
2.  **Linearization**: This substitution allows us to transform the nonlinear dispersion relation into a linear matrix eigenvalue problem of the form:
    ```math
    \mathcal{M} \mathbf{V} = \omega \mathbf{V}
    ```
    where $\mathcal{M}$ is a sparse matrix constructed from plasma parameters and auxiliary variables.
3.  **Global Solution**: By computing the eigenvalues of $\mathcal{M}$, we obtain **all** roots $\omega$ of the dispersion relation simultaneously.

**Key Advantages**:
-   **No Initial Guesses Required**: The solver finds all unstable and damped modes automatically.
-   **Robostness**: Effectively handles complex kinetic effects and arbitrary velocity distributions (via the [BO-HH](https://github.com/JuliaPlasma/PlasmaBO.jl) expansion).
-   **Efficiency**: Optimized for finding multiple modes in complex plasmas without manual intervention.

## Features

- Hermite-Hermite (BO-HH) expansion solver for arbitrary/analytic distributions
    - Maxwellian / BiMaxwellian
- Analytic Product-Bi-Kappa (BO-PBK) solver
- Multi-fluid solver
- Integration with [`ChargedParticles.jl`](https://juliaplasma.github.io/ChargedParticles.jl/dev/) for specifying different species, e.g., `Maxwellian("O-18 3+", n, Tpara)`

## Usage Examples

Check out the [ring beam instability example](ringbeam_Umeda12.md) for detailed usage instructions, also see [firehose instability example](firehose_Astfalk17.md) for using with arbitrary velocity distributions, [BO-PBK example](rlp_Cattaert07.md) for using with kappa distributions (BO-PBK), [cold plasma example](cold_plasma.md) for comparing kinetic and fluid solvers, and [dispersion surface tracking example](dispersion_surface.md) for 2D scanning and mode tracking.

### Solvers

```@docs; canonical = false
solve
```

BO-PBK ([`BOPBK`](@ref)) is an analytic, distribution-aware eigen-solver optimized for kappa plasmas, whereas BO-Arbitrary ([`BOHH`](@ref)) is a universal but numerically heavier framework that approximates any distribution at the cost of efficiency and low-κ accuracy.

```@docs; canonical = false
BOPBK
BOHH
BOFluid
```

## References

[xieBO20Plasma2021](@citet), [xiePDRFGeneralDispersion2014](@citet)

```@bibliography
```

## API Reference

```@index
```

```@autodocs
Modules = [PlasmaBO]
```
