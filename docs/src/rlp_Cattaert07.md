# Case: R-, L-, and P-mode waves

This page reproduces the “case 1” benchmark in BO-PBK (Cattaert et al., 2007) for the R-, L-, and P-mode branches.

The configuration is a single electron population with a strongly non-Maxwellian parallel kappa index (κ∥ = 1) and nearly Maxwellian perpendicular index (κ⟂ = 200), with oblique propagation at θ = 30°.

```@example rlp
using PlasmaBO
using PlasmaBO: q, me

B0 = 1.0e-6          # [Tesla]
θ = deg2rad(30)

n = 2.43e6           # [m^-3]
T = 2555.0           # [eV]

electron = Maxwellian(n, T; particle = :e)
param = HHSolverParam(electron, B0)

wce = abs(B0 * q / me)
ρce = param.vtp / wce
kn = 1 / ρce  # so that k/kn = kρce

kρ_scan = range(1.0e-4, 0.3; length = 80)
ks = kρ_scan .* kn

sol = solve_kinetic_dispersion([electron], B0, ks, θ; N =3)
```

```@example rlp
using CairoMakie

f, axs = plot(sol, kn, wce)
ylims!(axs[1], 0, 3)
f
```

