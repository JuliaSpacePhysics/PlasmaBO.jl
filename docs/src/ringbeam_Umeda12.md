# Case: Ring beam instability

This page demonstrates how to use the package to solve kinetic dispersion relations for the ring beam instability [umedaNumericalElectromagneticLinear2012](@citet).


### Simulation Parameters

- **Background Magnetic Field**: $B_0 = 96.24$ nT
- **Species**:
    1.  **Ring Beam Electrons** ($10\%$)
        -   Density: $n_r = 0.1 n_0$
        -   Temperature: $T_r = 51$ eV
        -   Drift Velocity: $v_{d,z} = 0.1 c$, $v_{d,\perp} = 0.05 c$
    2.  **Background Electrons** ($90\%$)
        -   Density: $n_b = 0.9 n_0$
        -   Temperature: $T_b = 51$ eV
        -   Drift Velocity: $v_{d} = 0$

```@example matrix
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed, debye_length
using Unitful

# Umeda 2012 ring beam configuration
B0 = 96.24e-9  # [Tesla]

e = Electron()
T = 51u"eV"
vth = thermal_speed(e, T)
plasma = Plasma(
    # Ring beam electrons (10% density): parallel drift `vd`, ring speed `vr`
    Population(e, GaussianRing(; vr_perp = 0.05u"c", w_perp = vth, vth_para = vth, vd_para = 0.1u"c"); n = 1e5),
    # Background electrons (90% density)
    Population(e, BiMaxwellian(; vth); n = 9e5);
    B0,
)

# Compute normalization
wce = abs(gyrofrequency(e, B0))
lambdaD = debye_length(e, 1e6, T)   # of the total electron density
kn = 1 / lambdaD

# Wave vector: k*λD = 0.03, θ = 40°
k = 0.03 / lambdaD
θ = deg2rad(40)
kx = k * sin(θ)
kz = k * cos(θ)

# J=12 provides good accuracy (J-pole approximation order)
ωs = solve(plasma, Wavevector(kx, kz), BOHH(; N=6, J=12))

# Filter for unstable modes (ω/ωce) with positive growth rate
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001*wce, ωs)[1] ./wce
```

## Dispersion Curve Scan

Scan k*λD from 0.01 to 0.3

```@repl matrix
k_ranges = (0.01:0.005:0.3) .* kn;
sol = solve(plasma, AngleSweep(; k = k_ranges, theta = θ), BOHH(; N=6))
```

```@example matrix
# k, ω pairs for initial branch points (see `BranchPoint` for more control over tracking)
initial_points = [
    (0.1 * kn, 0.3im * wce),
    (0.1 * kn, 0.1im * wce),
    (0.2 * kn, 0.25im * wce)
]

branches = track.(sol, initial_points)

# Extract individual branches
for (i, (k_branch, ω_branch)) in enumerate(branches)
    println("Branch $i:")
    println("  k range: $(minimum(k_branch)) to $(maximum(k_branch))")
    println("  Max growth rate: γ = $(maximum(imag.(ω_branch)))")
end
```

Plot all branches

```@example matrix
using CairoMakie

plot_branches(branches, kn, wce)
```