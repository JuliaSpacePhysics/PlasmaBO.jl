# Case: Dispersion surface tracking (2D scan)

This page demonstrates how to construct a dispersion *surface* by scanning a 2D parameter space (wave number `k` and propagation angle `θ`) and using branch tracking to follow a consistent mode across the scan.

We consider a typical electron–proton plasma and scan:

- `k` in normalized units `k/k_n`, where `k_n = ω_pi/c`
- `θ` from 5° to 90°

```@example surface
using PlasmaBO
using PlasmaPopulations: gyrofrequency, plasma_frequency
using PlasmaPopulations: C_SI as c0

B0 = 5e-9 # [Tesla]

n = 5.0e6
T = 12.94

p, e = Proton(), Electron()
plasma = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)

wn = abs(gyrofrequency(p, B0))
wpi = plasma_frequency(p, n)
kn = wpi / c0
wci = wn

# Kinetic solver settings (adjust upward for accuracy)
N = 3
```


## Dispersion Curve Scan

```@repl surface
using CairoMakie

θ = deg2rad(45);
sweep = AngleSweep(; k = collect(0.01:1:100) .* kn, theta = θ);
```

Every eigenvalue over the sweep, then a seed near a root at `k = 50 kn` picks the branch to follow:

```@example surface
survey = solve(DispersionProblem(plasma, sweep), BOHH(; N))
at = Wavevector(; k = 50kn, theta = θ)
branch = track(survey, (640 - 672im) * wn, at)
dispersion_diagram(branch; k_unit = kn, omega_unit = wn, xlabel = "k/kn")
```

## Dispersion Surface Scan (2D)

```@repl surface
ks = collect(0.1:10:100.0) .* kn;
θs = deg2rad.(10.0:10.0:90.0);
# Seed at a reference k and angle; the branch is linked over the (k, θ) grid
at = Wavevector(; k = 90.5kn, theta = deg2rad(20))
surface2d = solve(DispersionProblem(plasma, Seed((1000 - 3000im) * wn, at), AngleSweep(; k = ks, theta = θs)), BOHH(; N))
```

```@example surface
let fig = Figure(size = (900, 400))
    for (i, (part, label)) in enumerate(((real, "Re ω/ωci"), (imag, "Im ω/ωci")))
        ax = Axis3(fig[1, i]; xlabel = "k/kn", ylabel = "θ [°]", zlabel = label)
        surface!(ax, ks ./ kn, rad2deg.(θs), part.(surface2d.omega) ./ wn)
    end
    fig
end
```
