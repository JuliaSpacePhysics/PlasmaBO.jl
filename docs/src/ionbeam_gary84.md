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

## Parametric Variations (Test Cases)

The paper systematically varies specific parameters to isolate and study the growth rates of different electromagnetic instabilities.

- **Beam Density Variations:** The beam to proton density ratio ($n_b/n_p$) is tested at values of 0.005, 0.01, and 0.02.
- **Main Ion Beta Variations:** The main ion beta ($\beta_m$) is tested at values of 0.20, 1.0, and 5.0.
- **Beam Temperature Variations:** To represent different ion populations (like diffuse ions), the beam to main ion temperature ratio ($T_b/T_m$) is tested at 1, 10, 100, and 1000.
- **Beam Temperature Anisotropy:** To study the effect of perpendicular free energy, the ratio of perpendicular to parallel beam temperature is tested at values of 1 and 10.
- **Drift Velocity Variations:** The relative beam-main component drift speed ($v_0/v_A$) is swept from zero up to approximately 35 across multiple instability scans. Specific drift velocities examined in detail include $v_0 = 5 v_m, 10 v_m, 15 v_m, 20 v_m$, and $30 v_m$.
- **Oblique Propagation Case:** Oblique wave vector properties and higher order cyclotron resonances are primarily examined at a high drift speed of $v_0 = 30 v_m = 21.2 v_A$. The growth rate for the left oblique branch is also specifically detailed at a wavenumber of $k a_i = 0.44$.

## Example Setup in PlasmaBO.jl

The following example demonstrates how to set up the basic case with a drift velocity of $v_0 = 10 v_A$ using `PlasmaBO.jl`.

```@example gary84
using PlasmaBO
using PlasmaBO: c0, μ0, mp, q

# Global Parameters
vA_c = 1.0e-4
beta_m = 1.0

# Arbitrary magnetic field value (e.g., typical solar wind or magnetosphere values)
B0 = 1e-8 # 10 nT

# Derived parameters
vA = vA_c * c0

# Calculate density to match required vA
ne = B0^2 / (vA^2 * μ0 * mp)

nm = 0.99 * ne
nb = 0.01 * ne

Tm_eV = (beta_m * B0^2 / (2 * μ0) / nm) / q
Tb_eV = 10.0 * Tm_eV
Te_eV = 1.0 * Tm_eV

# Using an example case: v0 = 10 v_A
v0 = 10.0 * vA

# Zero current condition: n_e v_e = n_m v_m + n_b v_b
v_e = (nb / ne) * v0

# Define distributions
main_ion = Maxwellian(:p, nm, Tm_eV)
beam_ion = Maxwellian(:p, nb, Tb_eV; vdz=v0)
electron = Maxwellian(:e, ne, Te_eV; vdz=v_e)

species = [main_ion, beam_ion, electron]

# Compute normalizations
wcp = q * B0 / mp

# Define a wave vector range (k*vA/wcp)
k_ranges = (0.01:0.02:1.0) .* (wcp / vA)
θ = 0.0 # Parallel propagation

# Solve the dispersion relation
sol = solve(species, B0, k_ranges, θ; N=6)
```

### Dispersion Curve Scan

Here is a quick scan plotting the dispersion curve using `CairoMakie`.

```@example gary84
using CairoMakie

let
    # Plot using normalization (k_norm, w_norm)
    f, (ax1, ax2) = plot(sol, wcp / vA, wcp)
    
    ax1.ylabel = "Re(ω) / ωcp"
    ax2.ylabel = "Im(ω) / ωcp"
    ax2.xlabel = "k vA / ωcp"
    
    f
end
```
