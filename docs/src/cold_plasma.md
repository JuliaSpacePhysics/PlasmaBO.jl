# Case: Cold plasma (kinetic vs fluid)

This page compares the eigenmodes of the kinetic solver (`BOHH`) and the multi-fluid solver (`BOFluid`) for an electron–proton plasma.

The input parameters correspond to the following species table:

| species | n (m⁻³) | T (eV) |
|---:|---:|---:|
| ions (p) | 8.7e6 | 2.857e-3 |
| electrons (e) | 8.7e6 | 2.857e-3 |

```@example coldplasma
using PlasmaBO
using PlasmaPopulations: gyrofrequency, plasma_frequency
using PlasmaPopulations: C_SI as c0

B0 = 100e-9 # [Tesla]

n = 8.7e6
T = 2.857e-3

p, e = Proton(), Electron()
plasma = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)

wci = gyrofrequency(p, B0)
wpi = plasma_frequency(p, n)
kn = wpi / c0
```

## Dispersion Curve Scan

We scan `k` at fixed propagation angle `theta = 60°`.

```@example coldplasma
theta = deg2rad(60)
sweep = AngleSweep(; k = collect(0.01:2.0:100.0) .* kn, theta)
# The fluid closure is set on the solver, e.g. `BOFluid(; gamma_perp = 2.0, gamma_para = 3.0)`
kinetic = solve(DispersionProblem(plasma, sweep), BOHH())
fluid = solve(DispersionProblem(plasma, sweep), BOFluid())
```

```@example coldplasma
using CairoMakie

fig, (axr, axi) = dispersion_diagram(kinetic; k_unit = kn, omega_unit = wci, xlabel = "k/kn", label = "BOHH (kinetic)")
dispersion_diagram!((axr, axi), fluid; k_unit = kn, omega_unit = wci, color = :red, marker = :cross, label = "BOFluid")
axislegend(axi; position = :lb)
fig
```

The solvers agree on the fast electron-scale waves. Near the ion cyclotron frequency they differ:

```@example coldplasma
ylims!(axr, 0, 1.2)
ylims!(axi, -0.3, 0.02)
fig
```

The fluid ion cyclotron wave flattens at `ωci` undamped. The kinetic solver instead gives a fan of damped roots spreading from `ωci` linearly in `k`, its pole representation of ion cyclotron damping.
