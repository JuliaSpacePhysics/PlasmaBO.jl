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

Based on the parameter sets defined in Gary et al. (1984), here are complete `PlasmaBO.jl` configurations to demonstrate each instability mode.

### Right-Hand Resonant Ion Beam Instability

The following reproduces Gary et al. (1984) Figure 1, showing all four drift velocities in a 2×2 layout with $\omega_r/\omega_{cp}$ and $\gamma/\omega_{cp}$ plotted together on each panel. The wave number is normalized by the thermal proton gyroradius.
Note that $\frac{r_L}{d_i} = \frac{v_{ti}}{v_A} = \sqrt{\beta/2}$. For comparing with the original paper Figure 1, we normalize the wave number by $r_L$ instead of $d_i$. Also note that $a_m$ is used for $r_L$ in the original paper.

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
k_norm = 1 / sqrt(2) * ωcp / vm
k_ranges = (0.01:0.01:0.2) .* (ωcp / vm)

v0_vals   = [0.0, 10.0, 20.0, 30.0] .* vm
v0_labels = [L"v_0 = 0",
             L"v_0 = 10\,v_m",
             L"v_0 = 20\,v_m",
             L"v_0 = 30\,v_m"]

# v0=0: real seed targets the stable R-wave (Alfvén branch)
# v0>0: imaginary seed hunts the growing beam-resonant mode
initial_points = [
    (0.1 * ωcp / vm, 0.1 * ωcp),
    (0.1 * ωcp / vm, 0.08im * ωcp),
    (0.1 * ωcp / vm, 0.08im * ωcp),
    (0.1 * ωcp / vm, 0.08im * ωcp),
]

fig = Figure(size = (1000, 600), fontsize = 20)

for (idx, (v0, label)) in enumerate(zip(v0_vals, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col], title = label)

    col == 1 ? (ax.ylabel = L"$\omega\,/\,\omega_{cp}$") : hideydecorations!(ax, grid = false)
    row == 2 ? (ax.xlabel = L"$k\,r_{L,main}$") : hidexdecorations!(ax, grid = false)

    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    main_ion = Maxwellian(:p, nm, Tm_eV; vdz = v0m / c0)
    beam_ion = Maxwellian(:p, nb, Tb_eV; vdz = v0b / c0)
    electron = Maxwellian(:e, ne, Te_eV)
    species  = (main_ion, beam_ion, electron)

    sol = solve(species, B0, k_ranges, θ; N = 6, J = 12)
    # Tracking from given seeding points
    k_branch, ω_branch = track(sol, initial_points[idx])
    scatterlines!(ax, k_branch ./ k_norm, real.(ω_branch) ./ ωcp;
        linewidth = 2, color = :royalblue)
    scatterlines!(ax, k_branch ./ k_norm, imag.(ω_branch) ./ ωcp;
        linewidth = 2, color = :orangered)

    xlims!(ax, 0, 0.2)
    ylims!(ax, 0, 0.4)
end

Legend(fig[3, 1:2],
    [LineElement(color = :royalblue, linestyle = :dot, linewidth = 2),
     LineElement(color = :orangered, linestyle = :dot, linewidth = 2)],
    [L"$\omega_r / \omega_{cp}$", L"$\gamma / \omega_{cp}$"],
    orientation = :horizontal, framevisible = false)

fig
```

### Right-hand Nonresonant Ion Beam Instability

Requires either a high drift velocity or a high beam density.

```@example gary84_mode2
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8
vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
ωcp = q * B0 / mp

# Higher beam density ratio (e.g., 5%)
nm = 0.95 * ne
nb = 0.05 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
# thermal speed of main ions
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Higher drift speed, momentum-conserving
v0 = 20.0 * vA
v0m = -nb * v0 / (nm + nb)
v0b = v0m + v0
v_e = (nb / ne) * v0
θ = 0.0 # Parallel propagation

main_ion = Maxwellian(:p, nm, Tm_eV; vdz=v0m/c0)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0b/c0)
electron = Maxwellian(:e, ne, Te_eV; vdz=-v_e/c0)
species = (main_ion, beam_ion, electron)

k_ranges = (0.01:0.02:1.0) .* (ωcp / vm)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, 1 / sqrt(2) * ωcp / vm, ωcp)
ax1.ylabel = "Re(ω) / ωcp"
ax2.ylabel = "Im(ω) / ωcp"
ax2.xlabel = L"$k\, r_{L,main}$"
f
```

#### 3. Left-hand Resonant Ion Beam Instability
Dominant for very hot, diffuse beam distributions.
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
# thermal speed of main ions
vm = sqrt(Tm_eV * q / mp)
# Very large beam temperature
Tb_eV = 1000.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Momentum-conserving drifts
v0 = 10.0 * vA
v0m = -nb * v0 / (nm + nb)
v0b = v0m + v0
v_e = (nb / ne) * v0
θ = 0.0 # Parallel propagation

main_ion = Maxwellian(:p, nm, Tm_eV; vdz=v0m/c0)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0b/c0)
electron = Maxwellian(:e, ne, Te_eV; vdz=-v_e/c0)
species = (main_ion, beam_ion, electron)

k_ranges = (0.01:0.02:1.0) .* (ωcp / vm)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, 1 / sqrt(2) * ωcp / vm, ωcp)
ax1.ylabel = "Re(ω) / ωcp"
ax2.ylabel = "Im(ω) / ωcp"
ax2.xlabel = L"$k\, r_{L,main}$"
f
```

#### 4. Ion Cyclotron Anisotropy Instability
Driven by temperature anisotropy rather than drift velocity.
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
# thermal speed of main ions
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Zero drift speed
θ = 0.0 # Parallel propagation

# Temperature anisotropy T_perp / T_parallel = 4
Tparallel_b = 1.0 * Tb_eV
Tperp_b = 4.0 * Tb_eV

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tparallel_b, Tperp_b)
electron = Maxwellian(:e, ne, Te_eV)
species = (main_ion, beam_ion, electron)

k_ranges = (0.01:0.05:1.0) .* (ωcp / vm)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, 1 / sqrt(2) * ωcp / vm, ωcp)
ax1.ylabel = "Re(ω) / ωcp"
ax2.ylabel = "Im(ω) / ωcp"
ax2.xlabel = L"$k\, r_{L,main}$"
f
```

#### 5. Oblique Instabilities (Right & Left)
Require very high drift velocities and oblique propagation angles.
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
# thermal speed of main ions
vm = sqrt(Tm_eV * q / mp)
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Corresponds to 30 v_m in the paper (for beta_m = 1)
v0 = 21.2 * vA
v0m = -nb * v0 / (nm + nb)
v0b = v0m + v0
v_e = (nb / ne) * v0
# Oblique propagation angle
θ = π / 4

main_ion = Maxwellian(:p, nm, Tm_eV; vdz=v0m/c0)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0b/c0)
electron = Maxwellian(:e, ne, Te_eV; vdz=-v_e/c0)
species = (main_ion, beam_ion, electron)

# Scanning over a highly oblique angle
k_ranges = (0.01:0.02:2.0) .* (ωcp / vm)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, 1 / sqrt(2) * ωcp / vm, ωcp)
ax1.ylabel = "Re(ω) / ωcp"
ax2.ylabel = "Im(ω) / ωcp"
ax2.xlabel = L"$k\, r_{L,main}$"
f
```
