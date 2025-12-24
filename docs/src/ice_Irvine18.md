# Case: Ion cyclotron emission

Ion cyclotron emission (ICE) driven by a ring ion beam distribution in a magnetized fusion device

```@example ice
using PlasmaBO
using PlasmaBO: q, kb, ε0, me, c0, mp
using PlasmaBO: gyrofrequency, Alfven_speed

using DelimitedFiles

# Magnetic field (Tesla)
B0 = 2.1

# Wave vector angle between k and B0
θ = deg2rad(89.5)

# Species table taken from the original BO MATLAB input (docs/src/bo.in)
fpath = pkgdir(PlasmaBO, "test/ice_Irvine18.in")
tbl = readdlm(fpath, Float64; skipstart = 1)

species = map(eachrow(tbl)) do row
    q_s, m_s, n_s, Tz_s, Tp_s, vdz_s, vdr_s = row[1:7]
    Maxwellian(q_s, m_s, n_s, Tz_s, Tp_s; vdz = vdz_s, vdr = vdr_s)
end

ωn = gyrofrequency(B0, species[1])
vA = Alfven_speed(B0, species)
kn = ωn / vA

# Scan k/kn (bo_setup.m)
ks = (9.5:0.025:11.5) .* kn

# Solve using the same default accuracy parameters as the MATLAB setup
results = solve_kinetic_dispersion(species, B0, ks, θ; N = 12, J = 4)

# Extract the most unstable eigenvalue at each k
ω_maxγ = map(results.ωs) do ωs
    ωf = filter(isfinite, ωs)
    ωf[argmax(imag.(ωf))]
end
```

```@example ice
using CairoMakie

function plot_results(results, kn, wn)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "k [k_n]", ylabel = "ω [ω_ci]")
    ax2 = Axis(f[1, 2], xlabel = "k [k_n]", ylabel = "Im(ω) [ω_ci]")
    for (k, ωs) in  zip(results.ks, results.ωs)
        k_temp = fill(k / kn, length(ωs))
        scatter!(ax, k_temp, real.(ωs) ./ wn, color = :blue, markersize = 5)
        scatter!(ax2, k_temp, imag.(ωs) ./ wn, color = :red, markersize = 5)
    end
    return f, (ax, ax2)
end

let
    f, (ax, ax2) = plot_results(results, kn, ωn)
    ylims!(ax, [0, 15])
    ylims!(ax2, [-0.5, 0.5])
    f
end
```
