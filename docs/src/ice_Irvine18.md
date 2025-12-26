# Case: Ion cyclotron emission

Ion cyclotron emission (ICE) driven by a ring ion beam distribution in a magnetized fusion device

```@example ice
using PlasmaBO
using PlasmaBO: gyrofrequency, Alfven_speed

using DelimitedFiles

# Magnetic field (Tesla)
B0 = 2.1

# Wave vector angle between k and B0
θ = deg2rad(89.5)

# Species table taken from the original BO MATLAB input
fpath = pkgdir(PlasmaBO, "test/ice_Irvine18.in")
readdlm(fpath)
```

## Dispersion Curve Scan

```@example ice
tbl = readdlm(fpath, Float64; skipstart = 1)
species = map(eachrow(tbl)) do row
    Z, A, n_s, Tz_s, Tp_s, vdz_s, vdr_s = row[1:7]
    Maxwellian(n_s, Tz_s, Tp_s; vdz = vdz_s, vdr = vdr_s, Z , A)
end

ωn = gyrofrequency(B0, species[1])
vA = Alfven_speed(B0, species)
kn = ωn / vA

ks = (9.5:0.025:11.5) .* kn
results = solve(species, B0, ks, θ; N = 12, J = 4)
```

```@example ice
using CairoMakie

f, (ax, ax2) = plot(results, kn, ωn)
ylims!(ax, [0, 15])
ylims!(ax2, [-0.5, 0.5])
f
```
