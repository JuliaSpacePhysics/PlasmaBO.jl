# Case: Firehose instability (kappa distribution)

This page demonstrates test firehose instability using Hermite-Hermite (HH) expansion for bi-kappa distribution. 

Here we first generate numerical samples of the distribution function using `gen_fv2d`. Then we use `hermite_expansion` to compute the Hermite expansion coefficients (this would also work for arbitrary distributions).

```@example firehose
using PlasmaBO
using PlasmaBO: q, kb, ε0, me, c0, mp

B0 = 0.1  # [Tesla]
θ = deg2rad(45)
n = 5.0e19

me_mp = 1 / 1836 # [proton mass]
electron = Maxwellian(-1.0, me_mp, n, 496.683)

κz = 5.5
κx = 5.5
proton = BiKappa2(5.e19, κz, κx, 1986.734, 993.367)
data = gen_fv2d(proton)
```

```@example firehose
alm = hermite_expansion(data.fv, data.vz, data.vx, data.vtz, data.vtx).alm
proton_param = HHSolverParam(proton, B0; alm = alm)

kn = 31.0613
k = kn / 4
wci = proton_param.wc
species = (proton_param, electron)

ωs = solve_kinetic_dispersion(species, B0, k .* sincos(θ)...; N = 2, J = 24)
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
println("Unstable modes (ω/ωci): ", ω_unstable ./ wci)
```

## Dispersion Curve Scan

```@example firehose
k_ranges = (0.05:0.02:0.5) .* kn
results = solve_kinetic_dispersion(species, B0, k_ranges, θ; N = 2, J = 24);
```

```@example firehose
using CairoMakie

function plot_results(results, kn, wn)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "k [k_n]", ylabel = "ω [ω_ci]")
    ax2 = Axis(f[1, 2], xlabel = "k [k_n]", ylabel = "Im(ω) [ω_ci]")
    for (k, ωs) in zip(results.ks, results.ωs)
        k_temp = fill(k / kn, length(ωs))
        scatter!(ax, k_temp, real.(ωs) ./ wn, color = :blue, markersize = 5)
        scatter!(ax2, k_temp, imag.(ωs) ./ wn, color = :red, markersize = 5)
    end
    return f, (ax, ax2)
end
let
    f, (ax, ax2) = plot_results(results, kn, wci)
    ylims!(ax, [0, 1])
    ylims!(ax2, [-0.1, 0.1])
    f
end
```
