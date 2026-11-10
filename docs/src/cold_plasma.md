# Case: Cold plasma (kinetic vs fluid)

This page compares the eigenmodes computed by the kinetic solver (`solve`) and the multi-fluid solver (`solve` with `BOFluid`) for an electron–proton plasma.

The input parameters correspond to the following species table:

| species | n (m⁻³) | T (eV) |
|---:|---:|---:|
| ions (p) | 8.7e6 | 2.857e-3 |
| electrons (e) | 8.7e6 | 2.857e-3 |

```@example coldplasma
using PlasmaBO
using PlasmaPopulations: gyrofrequency, plasma_frequency
using PlasmaPopulations: C_SI as c0

B0 = 100e-9 # [Tesla]

n = 8.7e6
T = 2.857e-3

p, e = Proton(), Electron()
plasma = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)

wci = gyrofrequency(p, B0)
wpi = plasma_frequency(p, n)
kn = wpi / c0
```

## Dispersion Curve Scan

We scan `k` at fixed propagation angle `theta = 60°`.

```@repl coldplasma
theta = deg2rad(60)
kn_scan = 0.01:2.0:100.0;
k = kn_scan .* kn;
# Every mode at each wavevector
modes(alg) = map(k) do x
    solve(DispersionProblem(plasma, Wavevector(; k = x, theta)), alg).omega
end
kinetic = modes(BOHH());
fluid_ωs = modes(BOFluid());
# The fluid closure is set on the solver, e.g. `BOFluid(; gamma_perp = 2.0, gamma_para = 3.0)`
```

```@example coldplasma
using CairoMakie

let
    fig = Figure(size = (900, 360))
    ax1 = Axis(fig[1, 1], xlabel = "k [k_n]", ylabel = "Re(ω)/ω_ci")
    ax2 = Axis(fig[1, 2], xlabel = "k [k_n]", ylabel = "Im(ω)/ω_ci")

    for (k_norm, ωs_k, ωs_f) in zip(kn_scan, kinetic, fluid_ωs)
        xk = fill(k_norm, length(ωs_k))
        xf = fill(k_norm, length(ωs_f))

        scatter!(ax1, xk, real.(ωs_k) ./ wci; color=:transparent, strokecolor = (:blue, 0.35), strokewidth=2, marker = :circle)
        scatter!(ax1, xf .+ 0.05, real.(ωs_f) ./ wci; color=:red, marker = :cross)

        scatter!(ax2, xk, imag.(ωs_k) ./ wci; color=:transparent, strokecolor = (:blue, 0.35), strokewidth=2, marker = :circle)
        scatter!(ax2, xf .+ 0.05, imag.(ωs_f) ./ wci; color=:red, marker = :cross)
    end
    fig
end
```

Red crosses: Fluid solver results.

Blue circles: Kinetic solver results.

There is a slight difference at large k for the ion cyclotron wave, which is damped due to kinetic effect.