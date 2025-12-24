# Case: Firehose instability

```@example firehose
using PlasmaBO
using PlasmaBO: q, kb, ε0, me, c0, mp
using MAT: matopen

B0 = 0.1  # [Tesla]
θ = deg2rad(45)
n = 5.0e19

me_mp = 1 / 1836 # [proton mass]
electron = Maxwellian(-1.0, me_mp, n, 496.683)
fpath = pkgdir(PlasmaBO, "test/firehose_Astfalk17_fvceff1.mat")
proton_param = matopen(fpath) do file
    fvc = read(file, "fvc")
    HHSolverParam(q, mp, n, B0, fvc["vtz"], fvc["vtp"], 0.0, 0.0, fvc["alm"])
end

kn = 31.0613
k = kn / 4

wci = proton_param.wc

species = (proton_param, electron)

ωs = solve_kinetic_dispersion(species, B0, k .* sincos(θ)...; N = 2, J = 24)
ω_unstable = filter(ω -> isfinite(ω) && imag(ω) > 0.001 * wci, ωs)
println("Unstable modes (ω/ωci): ", ω_unstable ./ wci)
```

```@example firehose
using CairoMakie

k_ranges = (0.05:0.02:0.5) .* kn
results = solve_kinetic_dispersion((proton_param, electron), B0, k_ranges, θ; N = 2, J = 24)
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
