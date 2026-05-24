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
\det(\mathbf{D}(ω, \mathbf{k})) = \left|\mathbf{K}(ω, \mathbf{k})+\left(k k-k^2 \mathbf{I}\right) \frac{c^2}{ω^2}\right|=0
```

where $\mathbf{D}$ is the dispersion tensor, $ω$ is the complex frequency, $\mathbf{k}$ is the wave vector, $\mathbf{K}=\mathbf{I}+\mathbf{Q}$, $\mathbf{Q}=-\frac{\mathbf{\sigma}}{i ω \epsilon_0}$, and

```math
\mathbf{\sigma}=-i \sum_s \frac{q_s^2 n_s}{m_s} \sum_{n=-\infty}^{\infty} \int_{-\infty}^{\infty} \int_0^{\infty} \frac{2 \pi v_{\perp} d v_{\perp} d v_{\|}}{ω-n ω_{c s}-k_{\|} v_{\|}} \mathbf{\Pi}_s
```

with

```math
\mathbf{\Pi}_s=\left[\begin{array}{ccc}
A_s \frac{n^2 v_{\perp}}{\mu_s^2} J_n^2 & i A_s \frac{n v_{\perp}}{\mu_s} J_n J_n^{\prime} & B_s \frac{n v_{\perp}}{\mu_s} J_n^2 \\
-i A_s \frac{n v_{\perp}}{\mu_s} J_n J_n^{\prime} & A_s v_{\perp} J_n^{\prime 2} & -i B_s v_{\perp} J_n J_n^{\prime} \\
A_s \frac{n v_{\|}}{\mu_s} J_n^2 & i A_s v_{\|} J_n J_n^{\prime} & B_s v_{\|} J_n^2
\end{array}\right],
```

where $\mu_s=\frac{k_{\perp} v_{\perp}}{ω_{c s}}, A_s=\left(1-\frac{k_{\|} v_{\|}}{ω}\right) \frac{\partial f_{f_0}}{\partial v_{\perp}}+\frac{k_{\|} v_{\perp}}{ω} \frac{\partial f_{s 0}}{\partial v_{\|}}$, and $B_s=\frac{n ω_{c s} v_{\|}}{ω v_{\perp}} \frac{\partial f_{s 0}}{\partial v_{\perp}} +\left(1-\frac{n ω_{ω s}}{ω}\right) \frac{\partial f_{s_0}}{\partial v_{\|}}$.

This formulation is valid for non-relativistic, arbitrary gyrotropic distributions.

**Input**:
- **Plasma Parameters**: For each species $s$, density $n_s$, temperature $T_s$ (anisotropic $T_{\parallel}, T_{\perp}$), drift velocity $v_{d,s}$, charge $q_s$, and mass $m_s$.
    - **Note**: When arbitrary velocity distributions are provided, temperature are not needed as they are encoded in the distribution.
- **Coordinates & Fields**:
    - $z$: direction parallel to the background magnetic field ($B_0$).
    - $x$: one perpendicular direction (any perpendicular direction is equivalent).
    - Velocities: `vz` (parallel), `vx` (perpendicular), thermal speeds `vtz`, `vtx`.
- **Wave Parameters**: Wave vector $\mathbf{k} = (k_\perp, 0, k_\parallel)$, where $θ$ is the propagation angle between $k$ and $B_0$, and $k_∥ = k \cos θ$, $k_⊥ = k \sinθ$.

**Output**:
- **Complex Frequency**: $ω = ω_r + i\gamma$.
    - $ω_r$: Real frequency (oscillation).
    - $\gamma$: Growth rate (instability if $\gamma > 0$, damping if $\gamma < 0$).

## Methodology

Solving the kinetic dispersion relation is traditionally challenging because:
1.  **Transcendental Equations**: The dispersion tensor involves the Plasma Dispersion Function $Z(\zeta)$, making the equation $D(ω) = 0$ transcendental and highly nonlinear.
2.  **Root Finding Difficulty**: Traditional iterative methods (like Newton-Raphson) require good initial guesses. Without them, it's easy to miss modes or converge to the wrong root.
3.  **Multimode Nature**: Plasmas support multiple wave modes simultaneously (e.g., Alfvén, Whistlers, Langmuir), making it hard to track them all.

**PlasmaBO.jl** overcomes these challenges using a **Matrix Eigenvalue Method** ([xieEfficientFrameworkSolving2025](@citet), [xiePDRKGeneralKinetic2016](@citet)):

1.  **Rational Approximation**: We approximate the Plasma Dispersion Function $Z(\zeta)$ using a multi-pole (J-pole) expansion:

    ```math
    Z(\zeta) \approx \sum_{j=1}^J \frac{b_j}{\zeta - c_j}
    ```

2.  **Linearization**: This substitution allows us to transform the nonlinear dispersion relation into a linear matrix eigenvalue problem of the form:

    ```math
    \mathcal{M} \mathbf{V} = ω \mathbf{V}
    ```

    where $\mathcal{M}$ is a sparse matrix constructed from plasma parameters and auxiliary variables.
3.  **Global Solution**: By computing the eigenvalues of $\mathcal{M}$, we obtain **all** roots $ω$ of the dispersion relation simultaneously.

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

Check out the [ring beam instability example](ringbeam_Umeda12.md) for detailed usage instructions, also see [firehose instability example](firehose_Astfalk17.md) for using with arbitrary velocity distributions, [BO-PBK example](rlp_Cattaert07.md) for using with kappa distributions (BO-PBK), [cold plasma example](cold_plasma.md) for comparing kinetic and fluid solvers, [dispersion surface tracking example](dispersion_surface.md) for 2D scanning and mode tracking, and [wave polarization and handedness example](demo_polarization.md) for analyzing wave modes in cold plasma.

## Solvers

```@docs; canonical = false
solve
```

BO-PBK ([`BOPBK`](@ref)) is an analytic, distribution-aware eigen-solver optimized for kappa plasmas, whereas BO-Arbitrary ([`BOHH`](@ref)) is a universal but numerically heavier framework that approximates any distribution at the cost of efficiency and low-κ accuracy.

`BO-PBK` assumes **product-bi-kappa (PBK)** structure from the start:

```math
f(v_\parallel, v_\perp)=f_\parallel^{\kappa_\parallel}(v_\parallel) f_\perp^{\kappa_\perp}(v_\perp)
```

and uses exact modified plasma dispersion functions $Z_{\kappa}$ and their closed-form rational expansions. Loss-cone, drift, anisotropy are built in analytically.

`BOHH` expands the distribution as

```math
f_s(v_\parallel,v_\perp)=\sum_{l,m} a_{lm} \rho_l(v_\parallel) u_m(v_\perp)
```

using Hermite–Hermite (HH) bases. Hermite bases are Maxwellian-centered and accuracy degrades notably for **low κ (strong suprathermal tails)**.

```@docs; canonical = false
PlasmaBO.PBK.BOPBK
BOHH
BOFluid
```

### When should you use each?

Use **BO-PBK** if the plasma is reasonably described by PBK and you care about κ-dependent instability thresholds, loss-cone physics, and you want **speed + robustness**.

Use **BO-Arbitrary** if the distribution is numerical (e.g., from spacecraft data) and cannot be parameterized analytically (e.g., has plateaus, shoulders, or holes). You are willing to trade efficiency for **universality**.

## References

[xieBO20Plasma2021](@citet), [xiePDRFGeneralDispersion2014](@citet)

```@bibliography
```

## API Reference

```@index
```

```@autodocs
Modules = [PlasmaBO, PlasmaBO.PBK]
```
