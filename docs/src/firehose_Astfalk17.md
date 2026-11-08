# Case: Firehose instability (kappa distribution)

This page demonstrates test firehose instability using Hermite-Hermite (HH) expansion for bi-kappa distribution.
`BOHH` expands every non-Maxwellian model into a `HermiteVDF` about its thermal speeds and drift. Use `HermiteVDF(model; w_perp, w_para, N_para, N_perp)` to control the expansion, and `HermiteVDF(GridVDF(v_perp, v_para, f))` for tabulated distributions. The series is truncated, so low κ needs care: see [BO-PBK](rlp_Cattaert07.md) for integer κ∥.

### Simulation Parameters

- **Background Magnetic Field**: $B_0 = 0.1$ T
- **Species**:
    1.  **Electrons** (Maxwellian)
        -   Density: $n_e = 5 \times 10^{19}$ m$^{-3}$
        -   Temperature: $T_e = 496.683$ eV
    2.  **Protons** (Bi-Kappa)
        -   Density: $n_p = 5 \times 10^{19}$ m$^{-3}$
        -   Parallel Temperature: $T_{\parallel p} = 1986.734$ eV
        -   Perpendicular Temperature: $T_{\perp p} = 993.367$ eV
        -   Kappa Indices: $\kappa_{\parallel} = \kappa_{\perp} = 5.5$

```@example firehose
using PlasmaBO
using PlasmaPopulations: gyrofrequency

B0 = 0.1  # [Tesla]
θ = deg2rad(45)
n = 5.0e19

p, e = Proton(), Electron()
electron = BiMaxwellian(e; T = 496.683)
proton = ProductBiKappa(p; T_perp = 993.367, T_para = 1986.734, kappa_perp = 5.5, kappa_para = 5.5)
```

```@example firehose
plasma = Plasma(Population(p, proton; n), Population(e, electron; n); B0)

kn = 31.0613
k = kn / 4
wci = gyrofrequency(p, B0)

ωs = solve(plasma, Wavevector(k .* sincos(θ)...), BOHH(; N = 2, J = 24))
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
println("Unstable modes (ω/ωci): ", ω_unstable ./ wci)
```

## Dispersion Curve Scan

```@repl firehose
k_ranges = (0.05:0.02:0.5) .* kn;
sol = solve(plasma, AngleSweep(; k = k_ranges, theta = θ), BOHH(; N = 2, J = 24));
```

```@example firehose
using CairoMakie

let
    f, (ax, ax2) = plot(sol, kn, wci)
    ylims!(ax, [0, 1])
    ylims!(ax2, [-0.1, 0.1])
    f
end
```
