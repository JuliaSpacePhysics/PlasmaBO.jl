# Case: Ring beam instability

This page demonstrates how to use the package to solve kinetic dispersion relations for the ring beam instability [umedaNumericalElectromagneticLinear2012](@citet).

```@example matrix
using PlasmaBO
using PlasmaBO: q, me
using Unitful

# Umeda 2012 ring beam configuration
B0 = 96.24e-9  # [Tesla]

T = 51u"eV"
# Ring beam electrons (10% density)
# The first argument (optional) indicates particle type of the distribution (by default we use `proton`)
ring_beam = Maxwellian(:e, 1e5, T; vdz=0.1, vdr=0.05)
# Background electrons (90% density)
background = Maxwellian(:e, 9e5, T)

species = [ring_beam, background]

# Compute normalization
wce = abs(B0 * q / me)
lambdaD = Debye_length(species)

# Wave vector: k*λD = 0.03, θ = 40°
k = 0.03 / lambdaD
θ = deg2rad(40)
kx = k * sin(θ)
kz = k * cos(θ)

# J=12 provides good accuracy (J-pole approximation order)
ωs = solve_kinetic_dispersion(species, B0, kx, kz; N=6, J=12)

# Filter for unstable modes (ω/ωce) with positive growth rate
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001*wce, ωs)[1] ./wce
```

## Dispersion Curve Scan

```@example matrix
# Scan k*λD from 0.01 to 0.3
k_ranges = (0.01:0.004:0.3) ./ lambdaD
results = solve_kinetic_dispersion(species, B0, k_ranges, θ; N=6)
```

```@example matrix
using CairoMakie

# k, ω pairs for initial branch points (see `BranchPoint` for more control over tracking)
initial_points = [
    (0.1 / lambdaD, 0.3im * wce),
    (0.1 / lambdaD, 0.1im * wce),
    (0.2 / lambdaD, 0.25im * wce)
]

branches = track_dispersion_branches(results, initial_points)

# Extract individual branches
for (i, (k_branch, ω_branch)) in enumerate(branches)
    println("Branch $i:")
    println("  k range: $(minimum(k_branch)) to $(maximum(k_branch))")
    println("  Max growth rate: γ = $(maximum(imag.(ω_branch)))")
end

# Plot all branches
let xlabel = L"k λ_D", fig = Figure()
    ax1 = Axis(fig[1, 1]; xlabel, ylabel = "ωᵣ / ωₙ",)
    ax2 = Axis(fig[1, 2]; xlabel, ylabel = "γ / ωₙ")
    for (i, (k, ω)) in enumerate(branches)
        lines!(ax1, k * lambdaD, real.(ω) ./ wce, label = "Branch $i")
        lines!(ax2, k * lambdaD, imag.(ω) ./ wce)
    end
    fig
end
```