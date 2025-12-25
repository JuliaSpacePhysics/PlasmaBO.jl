# Case: Firehose instability (kappa distribution)

This page demonstrates test firehose instability using Hermite-Hermite (HH) expansion for bi-kappa distribution. 

Here we first generate numerical samples of the distribution function using `gen_fv2d`. Then we use `hermite_expansion` to compute the Hermite expansion coefficients (this would also work for arbitrary distributions).

```@example firehose
using PlasmaBO

B0 = 0.1  # [Tesla]
θ = deg2rad(45)
n = 5.0e19

electron = Maxwellian(:e, n, 496.683)

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

ωs = solve(species, B0, k .* sincos(θ)...; N = 2, J = 24)
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
println("Unstable modes (ω/ωci): ", ω_unstable ./ wci)
```

## Dispersion Curve Scan

```@repl firehose
k_ranges = (0.05:0.02:0.5) .* kn;
sol = solve(species, B0, k_ranges, θ; N = 2, J = 24);
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
