# Case: Dispersion surface tracking (2D scan)

This page demonstrates how to construct a dispersion *surface* by scanning a 2D parameter space (wave number `k` and propagation angle `θ`) and using branch tracking to follow a consistent mode across the scan.

We consider a typical electron–proton plasma and scan:

- `k` in normalized units `k/k_n`, where `k_n = ω_pi/c`
- `θ` from 5° to 90°

```@example surface
using PlasmaBO
using PlasmaBO: c0

B0 = 5e-9 # [Tesla]

n = 5.0e6
T = 12.94

ion = Maxwellian(:p, n, T)
electron = Maxwellian(:e, n, T)
species = (ion, electron)

wn = abs(B0 * ion.q / ion.m)
wpi = plasma_frequency(ion.q, n, ion.m)
kn = wpi / c0
wci = wn

# Kinetic solver settings (adjust upward for accuracy)
N = 3
```


## Dispersion Curve Scan

```@repl surface
using CairoMakie

θ = deg2rad(45);
k_ranges = (0.01:1:100) .* kn;
results = solve(species, B0, k_ranges, θ; N)
```

```@example surface
plot(results, kn, wn)
```

```@example surface
using PlasmaBO: plot_branches

# k, ω pairs for initial branch points (see `BranchPoint` for more control over tracking)
initial_point = (50 * kn, -600 * wn *im)
branch = track(results, initial_point)
f, (ax1, ax2) = plot_branches((branch,), kn, wn)
ylims!(ax2, -3500,500)
f
```

## Dispersion Surface Scan (2D)

```@repl surface
ks = (0.1:10:100.0) .* kn;
θs = deg2rad.(10.0:10.0:90.0);
res2d = solve(species, B0, ks, θs; N)
```

```@example surface
plot(res2d, kn, wn)
```

```@example surface
# Choose a reference angle and reference k for seeding the tracked mode
seed = (90.5 * kn, deg2rad(20), (1000 - 3000im) * wn)
branch = track(res2d, seed)
plot_branches((branch,), kn, wn)
```
