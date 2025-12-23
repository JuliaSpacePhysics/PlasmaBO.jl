# Examples

This page demonstrates how to use the package to solve kinetic dispersion relations for plasmas with arbitrary velocity distributions.

## Matrix Eigenvalue Solver

The matrix eigenvalue method finds all wave modes simultaneously by transforming the dispersion relation into a matrix eigenvalue problem using J-pole approximation for the plasma dispersion function.

This approach is more efficient to find multiple modes at once, and doesn't require initial guesses for the root finder.

Here we use the ring beam configuration from Umeda 2012 [umedaNumericalElectromagneticLinear2012](@citet).

```@example matrix
using PlasmaBO
using PlasmaBO: q, kb, ε0, me, c0

# Umeda 2012 ring beam configuration
B0 = 96.24e-9  # [Tesla]

me_mp = 1/1836 # [proton mass]
T = 51 # [eV]
# Ring beam electrons (10% density)
ring_beam = Maxwellian(-1.0, me_mp, 1e5, T; vdz=0.1, vdr=0.05)
# Background electrons (90% density)
background = Maxwellian(-1.0, me_mp, 9e5, T)

species = [ring_beam, background]

# Compute normalization
wce = abs(B0 * q / me)
lambdaD = Debye_length(species)

# Wave vector: k*λD = 0.03, θ = 40°
k = 0.03 / lambdaD
θ = deg2rad(40)
kx = k * sin(θ)
kz = k * cos(θ)

# Solve using matrix eigenvalue method
# J=12 provides good accuracy (J-pole approximation order)
ωs = solve_kinetic_dispersion(species, B0, kx, kz; N=6, J=12)

# Filter for unstable modes (ω/ωce) with positive growth rate
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001*wce, ωs)[1] ./wce
```

### Dispersion Curve Scan

```@example matrix
# Scan k*λD from 0.01 to 0.3
k_ranges = (0.01:0.0025:0.3) ./ lambdaD
results = solve_kinetic_dispersion(species, B0, k_ranges, θ; N=6)
```

```@example matrix
using CairoMakie

initial_points = [
    BranchPoint(0.1 / lambdaD, 0.3im * wce),
    BranchPoint(0.1 / lambdaD, 0.1im * wce),
    BranchPoint(0.2 / lambdaD, 0.25im * wce)
]

branches = track_dispersion_branches(results, initial_points)

# Extract individual branches
for (i, (k_branch, ω_branch)) in enumerate(branches)
    println("Branch $i:")
    println("  k range: $(minimum(k_branch)) to $(maximum(k_branch))")
    println("  Max growth rate: γ = $(maximum(imag.(ω_branch)))")
end

# Plot all branches
let xlabel = "k*λD", fig = Figure()
    ax1 = Axis(fig[1, 1]; xlabel, ylabel = "ωᵣ / ωₙ",)
    ax2 = Axis(fig[1, 2]; xlabel, ylabel = "γ / ωₙ")
    for (i, (k, ω)) in enumerate(branches)
        lines!(ax1, k * lambdaD, real.(ω) ./ wce, label = "Branch $i")
        lines!(ax2, k * lambdaD, imag.(ω) ./ wce)
    end

    fig
end
```