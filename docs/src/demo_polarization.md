# Case: Wave polarization and handedness

This page demonstrates how to analyze wave polarization, handedness, and the polarization ratio for various modes in a cold magnetized plasma. We solve the linear dispersion relation using the multi-fluid solver `BOFluid` and verify the calculated frequencies against the analytical cold-plasma dispersion relation.

## Physics Setup and Equations

We consider a cold, magnetized electron–proton plasma with the background magnetic field $\mathbf{B}_0$ aligned along the $+z$-direction.
The analytical cold-plasma dispersion relations for right-hand ($R$) and left-hand ($L$) circularly polarized waves propagating parallel to the magnetic field are:

```math
D_R(ω) = 1 - \frac{c^2 k^2}{ω^2} - \frac{ω_{pe}^2}{ω(ω - |ω_{ce}|)} - \frac{ω_{pi}^2}{ω(ω + ω_{ci})} = 0
```

```math
D_L(ω) = 1 - \frac{c^2 k^2}{ω^2} - \frac{ω_{pe}^2}{ω(ω + |ω_{ce}|)} - \frac{ω_{pi}^2}{ω(ω - ω_{ci})} = 0
```

where:
- $ω_{pe}, ω_{pi}$ are the electron and ion plasma frequencies, respectively.
- $ω_{ce}, ω_{ci}$ are the signed electron and ion cyclotron frequencies.

## Wave Polarization Analysis

In the code block below, we:
1. Setup the physical parameters ($B_0 = 10^{-4}$ T, $n = 10^{11}$ m$^{-3}$, $T = 0$ eV).
2. Construct the dispersion matrix using the fluid solver (`BOFluid`).
3. Compute the eigenvalues and eigenvectors using `eigen!`.
4. Identify and physically describe each of the 14 modes (e.g., Whistler waves, Ion Cyclotron waves, Langmuir waves, and uncoupled gyro-resonances).
5. Verify the numerical solutions against the analytical dispersion functions $D_R(ω)$ and $D_L(ω)$.

```@example polarization
using PlasmaBO
using PlasmaBO: c0, ε0, qe, me, mp
using LinearAlgebra
using Printf
using Markdown

# 1. Physics Setup
B0 = 1.0e-4          # Background magnetic field along +z [Tesla]
n = 1.0e11           # Plasma density [m^-3]
T = 0.0              # Temperature [eV] - cold limit for fluid

# Parallel propagation (along +z direction)
k = 1.0
kx = 0.0
kz = k

# Local helper functions to extract fields from fluid eigenvectors
electric_field(v) = (v[end-5], v[end-4], v[end-3])
magnetic_field(v) = (v[end-2], v[end-1], v[end])

# Define species
e_vdf = Maxwellian(:e, n, T)
species = [e_vdf, FluidSpecies(:p, n, T)]

# Construct the dispersion matrix using BOFluid solver
M = dispersion_matrix(species, B0, kx, kz, BOFluid)

# Solve the eigenvalue problem (eigen!)
decomp = eigen!(M)
ωs = decomp.values
V = decomp.vectors
```

We analyze each of the 14 modes by examining their frequency, electric field norm, polarization ratio, and handedness, comparing them against analytical cold-plasma theory:

```@example polarization
# Analytical Frequencies and Dispersion relation functions
wpe2 = n * qe^2 / (ε0 * me)
wpi2 = n * qe^2 / (ε0 * mp)
wce = qe * B0 / me     # Signed (negative for electrons)
wci = qe * B0 / mp     # Signed (positive for protons)
wpe = sqrt(wpe2)
wpi = sqrt(wpi2)

D_R(ω) = 1.0 - (c0*k/ω)^2 - wpe2 / (ω * (ω - abs(wce))) - wpi2 / (ω * (ω + wci))
D_L(ω) = 1.0 - (c0*k/ω)^2 - wpe2 / (ω * (ω + abs(wce))) - wpi2 / (ω * (ω - wci))

header = ["Index", "Frequency ω (rad/s)", "E-field Norm", "Handedness", "Pol Ratio", "Residual D", "Status", "Physical Description"]
rows = Any[]

# Sort all 14 solutions by |real(ω)|
sorted_indices = sort(1:length(ωs), by = i -> abs(real(ωs[i])))

for i in sorted_indices
    ω = ωs[i]
    v = V[:, i]
    
    # Extract components and analyze polarization
    Ex, Ey, Ez = electric_field(v)
    Bx, By, Bz = magnetic_field(v)
    E_norm = abs(Ex)^2 + abs(Ey)^2 + abs(Ez)^2
    hand = handedness(v, ω, kx, kz; threshold = 1e-4)
    
    # Residual and status defaults
    residual = NaN
    status = "N/A"
    desc = ""
    
    if abs(ω) < 1e1
        if abs(Bz) > 0.9
            desc = "spurious ∇·B null mode (δBz)"
        elseif abs(v[1]) > 0.9 && abs(v[5]) > 0.9
            desc = "Static Neutral Plasma Perturbation (physical)"
        elseif abs(v[1]) > 0.9
            desc = "spurious Gauss's Law violation (static δne)"
        elseif abs(v[5]) > 0.9
            desc = "spurious Gauss's Law violation (static δnp)"
        else
            desc = "Static Neutral Plasma Perturbation"
        end
    else
        # Calculate residual to verify if the computed eigenvalue satisfies 
        # the cold-plasma dispersion relation for its physical polarization handedness
        if real(ω) >= 0
            if hand == :R
                residual = abs(D_R(ω))
            elseif hand == :L
                residual = abs(D_L(ω))
            end
        else
            # For negative frequencies, R-mode satisfies D_L(ω) = 0 and L-mode satisfies D_R(ω) = 0
            if hand == :R
                residual = abs(D_L(ω))
            elseif hand == :L
                residual = abs(D_R(ω))
            end
        end
        
        if !isnan(residual)
            status = residual < 1e-4 ? "VERIFIED" : "DISCREPANT"
        end
        
        if E_norm < 1e-16
            residual = NaN
            status = "N/A"
            if abs(abs(ω) - abs(wci)) < 1e2
                desc = "Proton Cyclotron Gyration (Resonance)"
            else
                desc = "Uncoupled Species Gyration"
            end
        elseif hand == :linear
            desc = "Langmuir Wave " * (real(ω) > 0 ? "(+z)" : "(-z)")
        else
            # EM Waves
            pol = hand == :R ? "R-mode" : "L-mode"
            dir = real(ω) > 0 ? "+z" : "-z"
            if abs(real(ω)) > 2e8
                desc = "EM Light Wave ($(pol), $(dir))"
            elseif abs(real(ω)) > 1e6
                desc = "Whistler Wave ($(pol), $(dir))"
            else
                desc = "Ion Cyclotron Wave ($(pol), $(dir))"
            end
        end
    end
    
    # Format residual string
    res_str = isnan(residual) ? "N/A" : @sprintf("%.2e", residual)
    pol_ratio = polarization_ratio(v, kx, kz)
    # Use absolute value for display; handle non-finite values gracefully
    display_ratio = isfinite(pol_ratio) ? @sprintf("%.2e", abs(pol_ratio)) : "NaN"
    
    push!(rows, [
        string(i),
        "`" * @sprintf("%.3e + %.3eim", real(ω), imag(ω)) * "`",
        @sprintf("%.2e", E_norm),
        string(hand),
        display_ratio,
        res_str,
        status,
        desc
    ])
end

Markdown.MD(Markdown.Table([header, rows...], [:r, :l, :r, :c, :r, :r, :c, :l]))
```

### Summary of Characteristic Frequencies

We can also inspect the physical characteristics and cyclotron/plasma frequencies of the system:

```@example polarization
summary_str = """ # hide
**Characteristic Frequencies:** # hide
- Electron plasma frequency (\$ω_{pe}\$): $(@sprintf("%.3e", wpe)) rad/s # hide
- Ion cyclotron frequency (\$ω_{ci}\$): $(@sprintf("%.3e", wci)) rad/s # hide
- Electron cyclotron frequency (\$|ω_{ce}|\$): $(@sprintf("%.3e", abs(wce))) rad/s # hide
- High-frequency EM wave speed (\$c k\$): $(@sprintf("%.3e", c0 * k)) rad/s # hide

*Note: Handedness is defined relative to the background magnetic field \$\\mathbf{B}_0\$ (\$+\\mathrm{z}\$) assuming the wave convention \$e^{i(k_z z - ω t)}\$.* # hide
""" # hide
Markdown.parse(summary_str) # hide
```
