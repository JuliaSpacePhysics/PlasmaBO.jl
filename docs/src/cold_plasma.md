# Case: Cold plasma (kinetic vs fluid)

This page demonstrates a typical cold plasma configuration (electron–proton plasma) and compares the eigenmodes computed by the kinetic solver (`solve_kinetic_dispersion`) and the multi-fluid solver (`solve_fluid_dispersion`).

The input parameters correspond to the following species table:

| species | n (m⁻³) | T (eV) |
|---:|---:|---:|
| ions (p) | 8.7e6 | 2.857e-3 |
| electrons (e) | 8.7e6 | 2.857e-3 |

We scan `k` at fixed propagation angle `θ = 60°`.

```@example coldplasma
using PlasmaBO
using PlasmaBO: q, me, mp, c0

B0 = 100e-9 # [Tesla]
θ = deg2rad(60)

n = 8.7e6
T = 2.857e-3

ion = Maxwellian(:p, n, T)
electron = Maxwellian(:e, n, T)
species = (ion, electron)
# For more control over the fluid model, we can use `FluidSpecies`
# electron_f = FluidSpecies(:e, n, T; gamma_z = 1.0, gamma_p = 1.0)

wci = abs(B0 * q / mp)
wpi = plasma_frequency(q, n, mp)
kn = wpi / c0
```

## Dispersion Curve Scan

```@repl coldplasma
kn_scan = 0.01:2.0:100.0;
ks = kn_scan .* kn;
kinetic = solve_kinetic_dispersion(species, B0, ks, θ);
fluid_ωs = solve_fluid_dispersion(species, B0, ks, θ);
```

```@example coldplasma
using CairoMakie

let
    fig = Figure(size = (900, 360))
    ax1 = Axis(fig[1, 1], xlabel = "k [k_n]", ylabel = "Re(ω)/ω_ci")
    ax2 = Axis(fig[1, 2], xlabel = "k [k_n]", ylabel = "Im(ω)/ω_ci")

    for (k_norm, ωs_k, ωs_f) in zip(kn_scan, kinetic.ωs, fluid_ωs.ωs)
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