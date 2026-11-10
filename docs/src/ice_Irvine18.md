# Case: Ion cyclotron emission

Ion cyclotron emission (ICE) driven by a ring ion beam distribution

```@example ice
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed, alfven_speed
using PlasmaPopulations: E_SI as qe, MP_SI as mp, C_SI as c0

using DelimitedFiles

# Magnetic field (Tesla)
B0 = 2.1

# Species table taken from the original BO MATLAB input
fpath = pkgdir(PlasmaBO, "test/ice_Irvine18.in")
readdlm(fpath)
```

## Dispersion Curve Scan

```@example ice
tbl = readdlm(fpath, Float64; skipstart = 1)
# charge in e, mass in proton masses, temperatures in eV, drift and ring speeds in c
species = map(eachrow(tbl)) do row
    z, m, n, Tz, Tp, vd, vr = row[1:7]
    p = Particle(z * qe, m * mp)
    vdf = GaussianRing(; vr_perp = vr * c0, w_perp = thermal_speed(p, Tp), vth_para = thermal_speed(p, Tz), vd_para = vd * c0)
    Population(p, vdf; n)
end
plasma = Plasma(species...; B0)

ωn = gyrofrequency(species[1].particle, B0)
vA = alfven_speed(plasma)

theta = deg2rad(89.5)
kn = ωn / vA

k = (9.5:0.025:11.5) .* kn
box = ((0.0 - 0.5im) * ωn, (15.0 + 0.5im) * ωn)
results = solve(DispersionProblem(plasma, box, AngleSweep(; k = collect(k), theta)), BOHH(; N = 12, J = 4))
```

```@example ice
using CairoMakie

dispersion_diagram(results; k_unit = kn, omega_unit = ωn, xlabel = "k/kn")
```
