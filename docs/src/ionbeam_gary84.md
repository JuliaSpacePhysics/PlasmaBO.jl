# Ion Beam Instability (Gary et al. 1984)

This page reproduces the dispersion relations from the classical paper ([garyElectromagneticIonBeam1984](@citet)):
Gary, S. P., Smith, C. W., Lee, M. A., Goldstein, M. L., & Forslund, D. W. (1984). Electromagnetic ion beam instabilities. *Physics of Fluids*, 27(7), 1852-1862.

## Reference Distribution Function Parameters

These baseline parameters are used throughout the study unless otherwise specified.

**Global Parameters:**
- Main ion beta ($\beta_m$): 1.00
- Ratio of Alfvén speed to the speed of light ($v_A/c$): $1.0 \times 10^{-4}$

**Species Parameters:**

| Component | Density Ratio ($n_j/n_e$) | Mass Ratio ($m_j/m_e$) | Temp Ratio ($T_j/T_m$) | Anisotropy ($T_{\perp j}/T_{\parallel j}$) |
| :--- | :--- | :--- | :--- | :--- |
| Main ions (m) | 0.99 | 1836 | 1.00 | 1.00 |
| Beam ions (b) | 0.01 | 1836 | 10.00 | 1.00 |
| Electrons (e) | 1.00 | 1 | 1.00 | 1.00 |

### Summary of Electromagnetic Ion Beam Instabilities

| Instability Name | Propagation | Polarization | Resonance Type | Dominant Conditions |
| :--- | :--- | :--- | :--- | :--- |
| Right-hand resonant | Parallel | Right-hand circular (P = +1) | Beam resonant ($\|ζ_b\| \approx 1$) | Modest $v_0/v_A$, $T_b/T_m$, and $n_b/n_p$; resonant at large $k$ |
| Right-hand nonresonant | Antiparallel | Right-hand circular (P = -1 for $\omega_r < 0$) | Nonresonant ($\|ζ_j\| > 1$) | Sufficiently large $n_b/n_p$ or $v_0/v_A$; $T_b > T_m$ |
| Left-hand resonant | Parallel | Left-hand circular | Beam resonant ($\|ζ_b\| < 1$) | Large $T_b/T_m$ (diffuse distributions); lower threshold than right-hand nonresonant |
| Ion cyclotron anisotropy | Antiparallel | Left-hand circular | Resonant | Driven by temperature anisotropy ($T_{\perp}/T_{\parallel} > 1$, e.g., 4 to 9); modest $v_0/v_A$; persists at $v_0=0$ |
| Right oblique | Oblique to $\mathbf{B}$ | Right-hand elliptic (P > 0) | Higher harmonic resonant ($m = 2, 3...$) | High $v_0/v_A$ (e.g., $v_0 = 30 v_m$); free energy well above threshold of field-aligned mode |
| Left oblique | Oblique to $\mathbf{B}$ | Left-hand elliptic (P < 0) | Higher harmonic resonant ($m = 2, 3...$) | High drift speed ($v_0 = 30 v_m$); faster growth than right oblique $m=2$ |

## Simulating Specific Modes in PlasmaBO.jl

Based on the parameter sets defined in [garyElectromagneticIonBeam1984](@citet), here are complete `PlasmaBO.jl` configurations to demonstrate each instability mode.

### Right-Hand Resonant Ion Beam Instability

The following reproduces Figure 1, showing all four drift velocities in a 2×2 layout with $\omega_r/\omega_{cp}$ and $\gamma/\omega_{cp}$ plotted together on each panel. The wave number is normalized by the thermal proton gyroradius $r_L$ (i.e. $a_m$ in the original paper). Note that $\frac{r_L}{d_i} = \frac{v_{ti}}{v_A} = \sqrt{\beta/2}$.

```@example gary84_fig1
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8 # 10 nT

vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
vm = sqrt(Tm_eV * q / mp)  # main ion thermal speed
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

θ = 0.0
k_norm = ωcp / vm
# Note that the step size of k needs to be small enough
# to accurately capture the dispersion relation.
k_ranges = (0.01:0.001:0.2) .* (ωcp / vm)

v0_vals   = [0.0, 10.0, 20.0, 30.0] .* vm
v0_labels = [L"v_0 = 0",
             L"v_0 = 10\,v_m",
             L"v_0 = 20\,v_m",
             L"v_0 = 30\,v_m"]

# v0=0: real seed targets the stable R-wave (Alfvén branch)
# v0>0: imaginary seed hunts the growing beam-resonant mode
initial_points = [
    (0.1 * ωcp / vm, 0.15 * ωcp),
    (0.1 * ωcp / vm, 0.08im * ωcp),
    (0.05 * ωcp / vm, 0.1im * ωcp),
    (0.04 * ωcp / vm, 0.15im * ωcp),
]

sols = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    main_ion = Maxwellian(:p, nm, Tm_eV; vdz = v0m / c0)
    beam_ion = Maxwellian(:p, nb, Tb_eV; vdz = v0b / c0)
    electron = Maxwellian(:e, ne, Te_eV)
    species = (main_ion, beam_ion, electron)
    solve(species, B0, k_ranges, θ; N = 3)
end

fig = Figure(size = (900, 600), fontsize = 20)

for (idx, (sol, label)) in enumerate(zip(sols, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col], title = label)

    col == 1 ? (ax.ylabel = L"$\omega\,/\,\omega_{cp}$") : hideydecorations!(ax, grid = false)
    row == 2 ? (ax.xlabel = L"$k\,r_{L,main}$") : hidexdecorations!(ax, grid = false)

    # Tracking from given seeding points
    k_branch, ω_branch = track(sol, initial_points[idx])
    scatterlines!(ax, k_branch ./ k_norm, real.(ω_branch) ./ ωcp;
        linewidth = 2, color = :royalblue)
    scatterlines!(ax, k_branch ./ k_norm, imag.(ω_branch) ./ ωcp;
        linewidth = 2, color = :orangered)

    xlims!(ax, 0, 0.2)
    ylims!(ax, 0, 0.4)
end

colgap!(fig.layout, 30)

Legend(fig[3, 1:2],
    [LineElement(color = :royalblue, linestyle = :dot, linewidth = 2),
     LineElement(color = :orangered, linestyle = :dot, linewidth = 2)],
    [L"$\omega_r / \omega_{cp}$", L"$\gamma / \omega_{cp}$"],
    orientation = :horizontal, framevisible = false)

fig
```

### Right-hand Nonresonant Ion Beam Instability

Requires either a high drift velocity or a high beam density.

```@example gary84_fig2
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8 # 10 nT

vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
vm = sqrt(Tm_eV * q / mp)  # main ion thermal speed
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

θ = 0.0
k_norm = ωcp / vm
k_ranges = (0.01:0.001:0.2) .* k_norm

v0_vals   = [0.0, 10.0, 20.0, 30.0] .* vm
v0_labels = [L"v_0 = 0",
             L"v_0 = 10\,v_m",
             L"v_0 = 20\,v_m",
             L"v_0 = 30\,v_m"]

# The nonresonant mode is backward-propagating (ω_r < 0) and grows at low k.
# Seed with imaginary ω to target the growing nonresonant branch.
# For v0=0 there is no beam-driven nonresonant mode; target the stable Alfvén branch.
initial_points = [
    (0.1 * k_norm,  0.15 * ωcp),     # v0=0: stable Alfvén branch (real seed)
    (0.1 * k_norm, 0.1 * ωcp + 0.001im * ωcp),  # v0=10 vm: nonresonant mode (low k)
    (0.02 * k_norm, 0.0 + 0.1im * ωcp),   # v0=20 vm: stronger nonresonant growth
    (0.02 * k_norm, 0.0 + 0.25im * ωcp), # v0=30 vm: peak nonresonant instability
]

sols = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    main_ion = Maxwellian(:p, nm, Tm_eV; vdz = v0m / c0)
    beam_ion = Maxwellian(:p, nb, Tb_eV; vdz = v0b / c0)
    electron = Maxwellian(:e, ne, Te_eV)
    species = (main_ion, beam_ion, electron)
    solve(species, B0, k_ranges, θ; N = 3)
end

fig = Figure(size = (900, 600), fontsize = 20)

for (idx, (sol, label)) in enumerate(zip(sols, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col], title = label)

    col == 1 ? (ax.ylabel = L"$\omega\,/\,\omega_{cp}$") : hideydecorations!(ax, grid = false)
    row == 2 ? (ax.xlabel = L"$k\,r_{L,main}$") : hidexdecorations!(ax, grid = false)

    k_branch, ω_branch = track(sol, initial_points[idx])
    scatterlines!(ax, k_branch ./ k_norm, real.(ω_branch) ./ ωcp;
        linewidth = 2, color = :royalblue)
    scatterlines!(ax, k_branch ./ k_norm, imag.(ω_branch) ./ ωcp;
        linewidth = 2, color = :orangered)

    xlims!(ax, 0, 0.2)
    ylims!(ax, -0.2, 0.2)
end

colgap!(fig.layout, 30)

Legend(fig[3, 1:2],
    [LineElement(color = :royalblue, linestyle = :dot, linewidth = 2),
     LineElement(color = :orangered, linestyle = :dot, linewidth = 2)],
    [L"$\omega_r / \omega_{cp}$", L"$\gamma / \omega_{cp}$"],
    orientation = :horizontal, framevisible = false)

fig
```

#### 3. Left-hand Resonant Ion Beam Instability

This instability operates efficiently for hot, diffuse beam distributions where $T_b = 100 T_m$.
The setup generates the 4-panel layout corresponding to Figure 7.
A secondary axis is used because the growth rate $\gamma$ is scaled an order of magnitude smaller than the real frequency.

```@example gary84_mode3
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8

vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 100.0 * Tm_eV  # Diffuse beam distribution
Te_eV = 1.0 * Tm_eV

θ = 0.0
k_norm = ωcp / vm 
k_ranges = (0.01:0.002:0.25) .* k_norm

v0_vals   = [5.0, 10.0, 15.0, 20.0] .* vm
v0_labels = [L"v_0 = 5\,v_m", L"v_0 = 10\,v_m", L"v_0 = 15\,v_m", L"v_0 = 20\,v_m"]

initial_points = [
    (0.1 * k_norm, 0.1 * ωcp + 0.001im * ωcp),
    (0.1 * k_norm, 0.1 * ωcp + 0.005im * ωcp),
    (0.1 * k_norm, 0.1 * ωcp + 0.005im * ωcp),
    (0.1 * k_norm, 0.1 * ωcp + 0.005im * ωcp),
]

sols = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    main_ion = Maxwellian(:p, nm, Tm_eV; vdz = v0m / c0)
    beam_ion = Maxwellian(:p, nb, Tb_eV; vdz = v0b / c0)
    electron = Maxwellian(:e, ne, Te_eV)
    species = (main_ion, beam_ion, electron)
    solve(species, B0, k_ranges, θ; N = 4)
end

fig = Figure(size = (900, 600), fontsize = 18)

for (idx, (sol, label)) in enumerate(zip(sols, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax1 = Axis(fig[row, col], title = label)
    ax2 = Axis(fig[row, col], yaxisposition = :right)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    
    col == 1 ? (ax1.ylabel = L"$\omega_r\,/\,\omega_{cp}$") : hideydecorations!(ax1, grid = false)
    col == 2 ? (ax2.ylabel = L"$\gamma\,/\,\omega_{cp}$") : hideydecorations!(ax2, grid = false)
    row == 2 ? (ax1.xlabel = L"$k\,a_{m}$") : hidexdecorations!(ax1, grid = false)

    k_branch, ω_branch = track(sol, initial_points[idx])
    
    lines!(ax1, k_branch ./ k_norm, real.(ω_branch) ./ ωcp; linewidth = 2, color = :royalblue)
    scatter!(ax2, k_branch ./ k_norm, imag.(ω_branch) ./ ωcp; markersize = 6, color = :orangered)

    xlims!(ax1, 0, 0.2)
    xlims!(ax2, 0, 0.2)
    ylims!(ax1, 0, 0.2)
    ylims!(ax2, 0, 0.01) # Mode 3 growth rates are relatively small
end

colgap!(fig.layout, 30)
fig
```

#### 4. Ion Cyclotron Anisotropy Instability

Driven by temperature anisotropy rather than drift velocity. This mode persists even when the bulk drift velocity is zero.

```@example gary84_mode4
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8
vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

θ = 0.0

# Driving the instability via temperature anisotropy
Tparallel_b = 1.0 * Tb_eV
Tperp_b = 10.0 * Tb_eV

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tparallel_b, Tperp_b)
electron = Maxwellian(:e, ne, Te_eV)
species = (main_ion, beam_ion, electron)

k_norm = ωcp / vm 
k_ranges = (0.01:0.01:1.0) .* k_norm
sol = solve(species, B0, k_ranges, θ; N = 6)

fig = Figure(size = (800, 400), fontsize = 18)
ax1 = Axis(fig[1, 1], ylabel = L"$\omega_r / \omega_{cp}$", xlabel = L"$k\,a_{m}$")
ax2 = Axis(fig[1, 2], ylabel = L"$\gamma / \omega_{cp}$", xlabel = L"$k\,a_{m}$")

k_branch, ω_branch = track(sol, (0.2 * k_norm, 0.1im * ωcp))

lines!(ax1, k_branch ./ k_norm, real.(ω_branch) ./ ωcp, linewidth = 2, color = :royalblue)
lines!(ax2, k_branch ./ k_norm, imag.(ω_branch) ./ ωcp, linewidth = 2, color = :orangered)

fig
```

#### 5. Oblique Instabilities (Right & Left)

When examining oblique propagation at a large drift speed ($v_0 = 30 v_m$), higher-order cyclotron resonances emerge. This calculates the growth rate against the propagation angle $\theta$ for the Left Oblique mode at $k a_m = 0.44$, recreating Figure 10 which highlights the $m=2$ peak.

```@example gary84_mode5
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8

vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

v0 = 30.0 * vm
v0m = -nb * v0 / (nm + nb)
v0b = v0m + v0

main_ion = Maxwellian(:p, nm, Tm_eV; vdz = v0m / c0)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz = v0b / c0)
electron = Maxwellian(:e, ne, Te_eV)
species = (main_ion, beam_ion, electron)

k_target = 0.44 * ωcp / vm
θ_vals = range(66, 84, length = 40) .* (π / 180)

γ_vals = map(θ_vals) do θ
    sol = solve(species, B0, [k_target], θ; N = 6)
    # Extract the mode with the maximum growth rate at this configuration
    maximum(imag.(sol.ωs[1]))
end

fig = Figure(size = (600, 400), fontsize = 18)
ax = Axis(fig[1, 1], xlabel = L"$\theta$ (degrees)", ylabel = L"$\gamma / \omega_{cp}$", 
          title = "Left Oblique Mode", xminorticksvisible = true)

lines!(ax, θ_vals .* (180 / π), γ_vals ./ ωcp, linewidth = 2, color = :purple)

xlims!(ax, 65, 83)
ylims!(ax, 0.0, 0.08)
fig
```
