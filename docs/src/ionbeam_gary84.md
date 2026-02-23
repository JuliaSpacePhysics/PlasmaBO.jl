# Ion Beam Instability (Gary et al. 1984)

This page reproduces the dispersion relations from the classical paper:
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

#### 1. Right-hand Resonant Ion Beam Instability
This is the baseline case with modest beam velocity and temperature.
```@example gary84_mode1
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q
using CairoMakie

vA_c = 1.0e-4
beta_m = 1.0
B0 = 1e-8 # 10 nT

vA = vA_c * c0
ne = B0^2 / (vA^2 * μ0 * mp)
wcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Parameters targeting Right-hand Resonant mode
v0 = 10.0 * vA
v_e = (nb / ne) * v0
θ = 0.0 # Parallel propagation

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)
species = [main_ion, beam_ion, electron]

k_ranges = (0.01:0.02:1.0) .* (wcp / vA)
sol = solve(species, B0, k_ranges, θ;  N=6, J=12)

f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
ax1.ylabel = "Re(ω) / ωcp"; ax2.ylabel = "Im(ω) / ωcp"; ax2.xlabel = "k vA / ωcp"
f
```

#### 2. Right-hand Nonresonant Ion Beam Instability
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
wcp = q * B0 / mp

# Higher beam density ratio (e.g., 5%)
nm = 0.95 * ne
nb = 0.05 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Higher drift speed
v0 = 20.0 * vA
v_e = (nb / ne) * v0
θ = 0.0 # Parallel propagation

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)
species = [main_ion, beam_ion, electron]

k_ranges = (0.01:0.02:1.0) .* (wcp / vA)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
ax1.ylabel = "Re(ω) / ωcp"; ax2.ylabel = "Im(ω) / ωcp"; ax2.xlabel = "k vA / ωcp"
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
wcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
# Very large beam temperature
Tb_eV = 1000.0 * Tm_eV 
Te_eV = 1.0 * Tm_eV

v0 = 10.0 * vA
v_e = (nb / ne) * v0
θ = 0.0 # Parallel propagation

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)
species = [main_ion, beam_ion, electron]

k_ranges = (0.01:0.02:1.0) .* (wcp / vA)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
ax1.ylabel = "Re(ω) / ωcp"; ax2.ylabel = "Im(ω) / ωcp"; ax2.xlabel = "k vA / ωcp"
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
wcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Zero drift speed
v0 = 0.0 * vA
v_e = 0.0
θ = 0.0 # Parallel propagation

# Temperature anisotropy T_perp / T_parallel = 4
Tparallel_b = 1.0 * Tb_eV
Tperp_b = 4.0 * Tb_eV 

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tparallel_b, Tperp_b; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)
species = [main_ion, beam_ion, electron]

k_ranges = (0.01:0.02:1.0) .* (wcp / vA)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
ax1.ylabel = "Re(ω) / ωcp"; ax2.ylabel = "Im(ω) / ωcp"; ax2.xlabel = "k vA / ωcp"
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
wcp = q * B0 / mp

nm = 0.99 * ne
nb = 0.01 * ne
Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Corresponds to 30 v_m in the paper (for beta_m = 1)
v0 = 21.2 * vA 
v_e = (nb / ne) * v0
# Oblique propagation angle
θ = π / 4 

main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)
species = [main_ion, beam_ion, electron]

# Scanning over a highly oblique angle
k_ranges = (0.01:0.02:2.0) .* (wcp / vA)
sol = solve(species, B0, k_ranges, θ; N=6)

f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
ax1.ylabel = "Re(ω) / ωcp"; ax2.ylabel = "Im(ω) / ωcp"; ax2.xlabel = "k vA / ωcp"
f
```
