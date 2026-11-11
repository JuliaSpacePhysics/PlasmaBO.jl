# Case: Mirror mode

The benchmark is based on Gary (1993), p131, Fig.7.4, with

- **Propagation angle**: θ = 71°
- **Beta ratio**: $β_{i,∥}/β_{e,∥} = 1$
- **Temperature anisotropy**:
  - ions: $T_{i,⊥}/T_{i,∥} = 2$
  - electrons: $T_{e,⊥}/T_{e,∥} = 1$
- **Magnetic field**: $B_0 = 100e-9 T$

- **Species**:
  - Electron: n = 1e6, T∥ = 24840 eV, T⊥ = 24840 eV
  - Proton: n = 1e6, T∥ = 24840 eV, T⊥ = 49680 eV


```@example mirror
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed

B0 = 100e-9
theta = deg2rad(71)

n = 1e6
Tpara = 24840.0
Tperp = 49680.0

p, e = Proton(), Electron()
plasma = Plasma(
    Population(p, BiMaxwellian; T_perp = Tperp, T_para = Tpara, n),
    Population(e, BiMaxwellian; T = Tpara, n);
    B0,
)
ωₙ = gyrofrequency(p, B0)
ρᵢ = thermal_speed(p, Tperp) / sqrt(2) / ωₙ
k = (0.005:0.04:0.7) ./ ρᵢ

survey = solve(DispersionProblem(plasma, AngleSweep(; k = collect(k), theta)), BOHH())
# The most unstable root at each k
ωmax = argmax.(imag, survey.omega)
# Non-propagating: Re ω vanishes up to round-off
maximum(abs ∘ real, ωmax) / ωₙ
```

```@example mirror
using CairoMakie

scatterlines(k .* ρᵢ, imag.(ωmax) ./ ωₙ; axis = (; xlabel = "kρᵢ", ylabel = "γ / ωₙ"))
```