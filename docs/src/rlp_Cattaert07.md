# Case: R-, L-, and P-mode waves

This page reproduces the “case 1” benchmark ([cattaertObliquePropagationElectromagnetic2007](@citet)) in BO-PBK [baiBOPBKComprehensiveSolver2025](@citet) for the R-, L-, and P-mode branches.

The configuration is a single electron population with a strongly non-Maxwellian parallel kappa index (κ∥ = 1) and nearly Maxwellian perpendicular index (κ⟂ = 200), with oblique propagation at $θ = 30°$.

## Setup

```@example rlp
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed

B0 = 1.0e-6          # [Tesla]
θ = deg2rad(30)

n = 2.43e6           # [m^-3]
T = 2555.0           # [eV]
κz = 1.0
κx = 200.0

e = Electron()
vth = thermal_speed(e, T)
electron = Plasma(Population(e, ProductBiKappa(; vth, kappa_perp = κx, kappa_para = κz); n); B0)
wce = abs(gyrofrequency(e, B0))
ρce = sqrt(1 - 1 / κx) * vth / wce   # perpendicular kappa width θ⊥
kn = 1 / ρce  # so that k/kn = kρce

kρ_scan = range(1.0e-4, 0.3; length = 80)
ks = kρ_scan .* kn
```

## Analytic PBK (BO-PBK) solver

Solve the same case with the analytic PBK (BO-PBK) solver:

```@repl rlp
sweep = AngleSweep(; k = collect(ks), theta = θ);
box = ((0.0 - 1.0im) * wce, (3.0 + 7.0im) * wce);
sol_pbk = solve(DispersionProblem(electron, box, sweep), BOPBK())
```

```@example rlp
using CairoMakie

dispersion_diagram(sol_pbk; k_unit = kn, omega_unit = wce, xlabel = "kρce")
```

### Plot the tracked PBK branches:

Seeds at $kρ_{ce} = 0.1$ are used to identify the branches: upper R-X, O(P), L-X, lower R-X

```@example rlp
seed_kρ = 0.1
seeds_ω = [1.37366, 1.07271, 0.75519, 0.39861] .* wce

at = Wavevector(; k = seed_kρ * kn, theta = θ)
# Every eigenvalue, so the branches are not cut where they leave `box`
all_pbk = solve(DispersionProblem(electron, sweep), BOPBK())
branches_pbk = [track(all_pbk, ω0, at) for ω0 in seeds_ω]

f, (ax1, ax2) = dispersion_diagram(branches_pbk; k_unit = kn, omega_unit = wce, xlabel = "kρce")
ylims!(ax2, -0.01, 0.001)
f
```

## General HH solver

Solve the same case with the general Hermite-Hermite (HH) basis solver:

```@example rlp
using CairoMakie

sol_HH = solve(DispersionProblem(electron, box, sweep), BOHH(; N = 4, J = 16))
all_HH = solve(DispersionProblem(electron, sweep), BOHH(; N = 4, J = 16))
branches_HH = [track(all_HH, ω0, at) for ω0 in seeds_ω]

dispersion_diagram(sol_HH; k_unit = kn, omega_unit = wce, xlabel = "kρce")
```

Plot the tracked branches

```@example rlp
dispersion_diagram(branches_HH; k_unit = kn, omega_unit = wce, xlabel = "kρce")
```

We can see that the results from the HH solver is very different from the PBK solver. It does not yield the same branches as the PBK solver due to the small parallel kappa index.

As a comparison, we solve it for a Maxwellian electron population:

```@example rlp
maxwellian = Plasma(Population(e, BiMaxwellian(; vth); n); B0)
sol_M = solve(DispersionProblem(maxwellian, box, sweep), BOHH(; N = 3))

dispersion_diagram(sol_M; k_unit = kn, omega_unit = wce, xlabel = "kρce")
```