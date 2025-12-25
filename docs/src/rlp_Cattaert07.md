# Case: R-, L-, and P-mode waves

This page reproduces the “case 1” benchmark ([cattaertObliquePropagationElectromagnetic2007](@citet)) in BO-PBK [baiBOPBKComprehensiveSolver2025](@citet) for the R-, L-, and P-mode branches.

The configuration is a single electron population with a strongly non-Maxwellian parallel kappa index (κ∥ = 1) and nearly Maxwellian perpendicular index (κ⟂ = 200), with oblique propagation at $θ = 30°$.

## Setup

```@example rlp
using PlasmaBO
using PlasmaBO: q, me

B0 = 1.0e-6          # [Tesla]
θ = deg2rad(30)

n = 2.43e6           # [m^-3]
T = 2555.0           # [eV]
κz = 1.0
κx = 200.0

electron = BiKappa2(:e, n, κz, κx, T; sigma = 0.0)
wce = abs(B0 * q / me)
ρce = electron.vtp / wce
kn = 1 / ρce  # so that k/kn = kρce

kρ_scan = range(1.0e-4, 0.3; length = 80)
ks = kρ_scan .* kn
```

## Analytic PBK (BO-PBK) solver

Solve the same case with the analytic PBK (BO-PBK) solver:

```@repl rlp
sol_pbk = solve(electron, B0, ks, θ, BOPBK);
```

```@example rlp
using CairoMakie

f, axs = plot(sol_pbk, kn, wce)
ylims!(axs[1], 0, 3)
f
```

### Plot the tracked PBK branches:

Seeds at kρce = 0.1 are used to identify the branches: upper R-X, O(P), L-X, lower R-X

```@example rlp
seed_kρ = 0.1
seeds_ω = [1.37366, 1.07271, 0.75519, 0.39861] .* wce

initial_points = [(seed_kρ * kn, ω0) for ω0 in seeds_ω]
branches_pbk = track.(sol_pbk, initial_points)

f, (ax1, ax2) = plot_branches(branches_pbk, kn, wce)
ylims!(ax1, nothing, 3)
ylims!(ax2, -0.01, 0.001)
f
```

## General HH solver

Solve the same case with the general Hermite-Hermite (HH) basis solver:

```@example rlp
using CairoMakie

sol_HH = solve(electron, B0, ks, θ; N = 4, J = 16)
branches_HH = track.(sol_HH, initial_points)

f, axs = plot(sol_HH, kn, wce)
ylims!(axs[1], 0, 3)
f
```

Plot the tracked branches

```@example rlp
plot_branches(branches_HH, kn, wce)
```

We can see that the results from the HH solver is very different than the PBK solver. It does not yield the same branches as the PBK solver due to the small parallel kappa index.

As a comparison, we solve it for a Maxwellian electron population:

```@example rlp
electron = Maxwellian(:e, n, T)
sol_M = solve([electron], B0, ks, θ; N = 3)

f, axs = plot(sol_M, kn, wce)
ylims!(axs[1], 0, 3)
f
```