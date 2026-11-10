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
| Right-hand resonant | With beam | Right-hand circular (P = +1) | Beam resonant ($\|ζ_b^+\| \lesssim 1$) | Largest $\gamma_m$ for modest $v_0/v_A$, $T_b/T_m$ and $n_b/n_p$ |
| Right-hand nonresonant | Against beam ($\omega_r < 0$) | Right-hand circular (P = +1) | Nonresonant ($\|ζ_j^-\| > 1$ for all $j$) | Largest $\gamma_m$ for sufficiently large $n_b/n_p$ or $v_0/v_A$ |
| Left-hand resonant | With beam | Left-hand circular (P = -1) | Beam resonant ($\|ζ_b^-\| < 1$) | Hot beams, beam thermal speed $\gtrsim v_{0b}$; much lower threshold than right-hand nonresonant |
| Ion cyclotron anisotropy | Against beam | Left-hand circular | Resonant | $T_{\perp b} > T_{\parallel b}$ (observed 4 to 9, Fig. 6 uses 10) and $v_0/v_A$ not too large; persists at $v_0 = 0$ |
| Right oblique | Oblique to $\mathbf{B}$ | Right-hand elliptic (P > 0) | Higher harmonic resonant ($m = 2, 3...$) | High $v_0/v_A$ (e.g., $v_0 = 30 v_m$); free energy well above the right-hand nonresonant threshold |
| Left oblique | Oblique to $\mathbf{B}$ | Evolves from the left-hand $\omega_r > 0$ wave; nearly linear at large $\theta$ | Higher harmonic resonant ($m = 2, 3...$) | High drift speed ($v_0 = 30 v_m$); its $m = 2$ peak grows faster than the right oblique one |

## Simulating Specific Modes in PlasmaBO.jl

Based on the parameter sets defined in [garyElectromagneticIonBeam1984](@citet), here are complete `PlasmaBO.jl` configurations to demonstrate each instability mode.

### Right-hand resonant ion beam instability (Fig. 1)

The following reproduces Figure 1, showing all four drift velocities in a 2×2 layout with $\omega_r/\omega_{cp}$ and $\gamma/\omega_{cp}$ plotted together on each panel. The wave number is normalized by the main-ion thermal gyroradius $a_m = v_m/\Omega_p$, with $v_m = \sqrt{T_m/m_p}$; $v_m/v_A = \sqrt{\beta_m/2}$.

```@example gary84_fig1
using PlasmaBO
using PlasmaPopulations: C_SI as c0, MU0_SI as μ0, MP_SI as mp, E_SI as q
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

# (k, ω) seeds: v0=0 on the stable R-wave (Alfvén branch), v0>0 on the growing beam-resonant mode
initial_points = [
    (0.1 * ωcp / vm, 0.15 * ωcp),
    (0.1 * ωcp / vm, (0.167 + 0.0855im) * ωcp),
    (0.05 * ωcp / vm, (0.127 + 0.14im) * ωcp),
    (0.04 * ωcp / vm, (0.187 + 0.148im) * ωcp),
]

plasmas = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    plasma = Plasma(
        Population(Proton(), BiMaxwellian; T = Tm_eV, vd_para = v0m, n = nm),
        Population(Proton(), BiMaxwellian; T = Tb_eV, vd_para = v0b, n = nb),
        Population(Electron(), BiMaxwellian; T = Te_eV, n = ne);
        B0,
    )
    plasma
end
alg = BOHH(; N = 3)

fig = Figure(size = (900, 600), fontsize = 20)

for (idx, (plasma, label)) in enumerate(zip(plasmas, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col], title = label)

    col == 1 ? (ax.ylabel = L"$\omega\,/\,\omega_{cp}$") : hideydecorations!(ax, grid = false)
    row == 2 ? (ax.xlabel = L"$k\,a_m$") : hidexdecorations!(ax, grid = false)

    # Tracking from given seeding points
    k0, ω0 = initial_points[idx]
    seed = Seed(ω0, Wavevector(; k = k0, theta = θ))
    ω_branch = solve(DispersionProblem(plasma, seed, AngleSweep(; k = collect(k_ranges), theta = θ)), alg).omega
    scatterlines!(ax, k_ranges ./ k_norm, real.(ω_branch) ./ ωcp;
        linewidth = 2, color = :royalblue)
    scatterlines!(ax, k_ranges ./ k_norm, imag.(ω_branch) ./ ωcp;
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

### Right-hand nonresonant ion beam instability (Fig. 2)

Requires either a high drift velocity or a high beam density. Figure 2 shows the two waves of Eq. (5−): the left-hand Alfvén wave ($\omega_r > 0$) and the right-hand wave propagating against the beam ($\omega_r < 0$). Above threshold they merge at small $k$ into a pair with equal $\omega_r$ and opposite $\gamma$; the growing one is the nonresonant instability.

```@example gary84_fig2
using PlasmaBO
using PlasmaPopulations: C_SI as c0, MU0_SI as μ0, MP_SI as mp, E_SI as q
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

# Seeds at k a_m = 0.18 on the two Eq. (5−) waves: ω_r > 0 left-hand, ω_r < 0 right-hand.
# The right-hand one carries the nonresonant instability at small k for v0 ≥ 20 vm.
k_seed = 0.18 * k_norm
seeds = [(0.1981, -0.3004), (0.1665, -0.2829), (0.1183, -0.2525), (0.0574, -0.2091)]

plasmas = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    plasma = Plasma(
        Population(Proton(), BiMaxwellian; T = Tm_eV, vd_para = v0m, n = nm),
        Population(Proton(), BiMaxwellian; T = Tb_eV, vd_para = v0b, n = nb),
        Population(Electron(), BiMaxwellian; T = Te_eV, n = ne);
        B0,
    )
    plasma
end
alg = BOHH(; N = 3)

fig = Figure(size = (900, 600), fontsize = 20)

for (idx, (plasma, label)) in enumerate(zip(plasmas, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col], title = label)

    col == 1 ? (ax.ylabel = L"$\omega\,/\,\omega_{cp}$") : hideydecorations!(ax, grid = false)
    row == 2 ? (ax.xlabel = L"$k\,a_m$") : hidexdecorations!(ax, grid = false)

    survey = solve(DispersionProblem(plasma, AngleSweep(; k = collect(k_ranges), theta = θ)), alg)
    for ω0 in seeds[idx]
        ω_branch = track(survey, ω0 * ωcp, Wavevector(; k = k_seed, theta = θ)).omega
        scatterlines!(ax, k_ranges ./ k_norm, real.(ω_branch) ./ ωcp;
            linewidth = 2, color = :royalblue)
        scatterlines!(ax, k_ranges ./ k_norm, imag.(ω_branch) ./ ωcp;
            linewidth = 2, color = :orangered)
    end

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

### Left-hand resonant ion beam instability (Fig. 7)

This instability operates efficiently for hot, diffuse beam distributions where $T_b = 100 T_m$.
The setup generates the 4-panel layout corresponding to Figure 7.
A secondary axis is used because the growth rate $\gamma$ is scaled an order of magnitude smaller than the real frequency.

```@example gary84_mode3
using PlasmaBO
using PlasmaPopulations: C_SI as c0, MU0_SI as μ0, MP_SI as mp, E_SI as q
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

plasmas = map(v0_vals) do v0
    v0m = -nb * v0 / (nm + nb)
    v0b = v0m + v0

    plasma = Plasma(
        Population(Proton(), BiMaxwellian; T = Tm_eV, vd_para = v0m, n = nm),
        Population(Proton(), BiMaxwellian; T = Tb_eV, vd_para = v0b, n = nb),
        Population(Electron(), BiMaxwellian; T = Te_eV, n = ne);
        B0,
    )
    plasma
end
alg = BOHH(; N = 4)

fig = Figure(size = (900, 600), fontsize = 18)

for (idx, (plasma, label)) in enumerate(zip(plasmas, v0_labels))
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax1 = Axis(fig[row, col], title = label)
    ax2 = Axis(fig[row, col], yaxisposition = :right)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    
    col == 1 ? (ax1.ylabel = L"$\omega_r\,/\,\omega_{cp}$") : hideydecorations!(ax1, grid = false)
    col == 2 ? (ax2.ylabel = L"$\gamma\,/\,\omega_{cp}$") : hideydecorations!(ax2, grid = false)
    row == 2 ? (ax1.xlabel = L"$k\,a_{m}$") : hidexdecorations!(ax1, grid = false)

    k0, ω0 = initial_points[idx]
    seed = Seed(ω0, Wavevector(; k = k0, theta = θ))
    ω_branch = solve(DispersionProblem(plasma, seed, AngleSweep(; k = collect(k_ranges), theta = θ)), alg).omega
    
    lines!(ax1, k_ranges ./ k_norm, real.(ω_branch) ./ ωcp; linewidth = 2, color = :royalblue)
    scatter!(ax2, k_ranges ./ k_norm, imag.(ω_branch) ./ ωcp; markersize = 6, color = :orangered)

    xlims!(ax1, 0, 0.2)
    xlims!(ax2, 0, 0.2)
    ylims!(ax1, 0, 0.2)
    ylims!(ax2, 0, 0.01) # Mode 3 growth rates are relatively small
end

colgap!(fig.layout, 30)
fig
```

### Ion cyclotron anisotropy instability (Fig. 6)

Driven by temperature anisotropy rather than drift velocity. This mode persists even when the bulk drift velocity is zero: with $T_{\perp b} = 10\,T_{\parallel b}$ and $v_0 = 0$, its maximum growth rate is the $v_0 = 0$ end of the hatched line in Fig. 6, $\gamma_m \approx 0.07\,\Omega_p$.

```@example gary84_mode4
using PlasmaBO
using PlasmaPopulations: C_SI as c0, MU0_SI as μ0, MP_SI as mp, E_SI as q
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

p = Proton()
plasma = Plasma(
    Population(p, BiMaxwellian; T = Tm_eV, n = nm),
    Population(p, BiMaxwellian; T_perp = Tperp_b, T_para = Tparallel_b, n = nb),
    Population(Electron(), BiMaxwellian; T = Te_eV, n = ne);
    B0,
)

k_norm = ωcp / vm 
k_ranges = (0.01:0.005:0.4) .* k_norm

fig = Figure(size = (800, 400), fontsize = 18)
ax1 = Axis(fig[1, 1], ylabel = L"$\omega_r / \omega_{cp}$", xlabel = L"$k\,a_{m}$")
ax2 = Axis(fig[1, 2], ylabel = L"$\gamma / \omega_{cp}$", xlabel = L"$k\,a_{m}$")

seed = Seed((-0.22 + 0.0661im) * ωcp, Wavevector(; k = 0.2k_norm, theta = θ))
ω_branch = solve(DispersionProblem(plasma, seed, AngleSweep(; k = collect(k_ranges), theta = θ)), BOHH(; N = 6)).omega

lines!(ax1, k_ranges ./ k_norm, real.(ω_branch) ./ ωcp, linewidth = 2, color = :royalblue)
lines!(ax2, k_ranges ./ k_norm, imag.(ω_branch) ./ ωcp, linewidth = 2, color = :orangered)

fig
```

### Oblique instabilities (Fig. 10)

When examining oblique propagation at a large drift speed ($v_0 = 30 v_m$), higher-order cyclotron resonances emerge. This calculates the growth rate against the propagation angle $\theta$ for the left oblique mode at $k a_m = 0.44$, recreating Figure 10: relative maxima at the $m = 3$ (near 76°) and $m = 2$ (near 81°) beam cyclotron resonances, $\omega_r \approx k_z v_{0b} - m\Omega_p$.

```@example gary84_mode5
using PlasmaBO
using PlasmaPopulations: C_SI as c0, MU0_SI as μ0, MP_SI as mp, E_SI as q
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

plasma = Plasma(
    Population(Proton(), BiMaxwellian; T = Tm_eV, vd_para = v0m, n = nm),
    Population(Proton(), BiMaxwellian; T = Tb_eV, vd_para = v0b, n = nb),
    Population(Electron(), BiMaxwellian; T = Te_eV, n = ne);
    B0,
)

k_target = 0.44 * ωcp / vm
θ_vals = range(66, 84, length = 40) .* (π / 180)

γ_vals = map(θ_vals) do θ
    sol = solve(DispersionProblem(plasma, Wavevector(; k = k_target, theta = θ)), BOHH(; N = 6))
    # The left oblique mode is the only growing one here
    maximum(imag, sol.omega)
end

fig = Figure(size = (600, 400), fontsize = 18)
ax = Axis(fig[1, 1], xlabel = L"$\theta$ (degrees)", ylabel = L"$\gamma / \omega_{cp}$", 
          title = "Left Oblique Mode", xminorticksvisible = true)

lines!(ax, θ_vals .* (180 / π), γ_vals ./ ωcp, linewidth = 2, color = :purple)

xlims!(ax, 65, 83)
ylims!(ax, 0.0, 0.08)
fig
```
