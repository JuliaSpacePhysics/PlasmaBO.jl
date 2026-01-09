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

B0 = 100e-9
θ = deg2rad(71)

n = 1e6
Tpara = 24840.0
Tperp = 49680.0

electron = Maxwellian(:e, n, Tpara)
proton   = Maxwellian(n, Tpara, Tperp)

species = (proton, electron)
params = HHSolverParam.(species, B0)
ρᵢ = params[1].ρc
ωₙ = params[1].wc
ks = (0.005:0.04:0.7) ./ ρᵢ

sol = solve(species, B0, ks, θ)

# Extract the most unstable mode at each k (by growth rate)
ωmax = vec(argmax.(imag, sol.ωs))
```

```@example mirror
using CairoMakie

let fig = Figure()
    ax1 = Axis(fig[1,1]; xlabel = "k*λD", ylabel = "ωᵣ / ωₙ")
    ax2 = Axis(fig[1,2]; xlabel = "k*λD", ylabel = "γ / ωₙ")
    scatterlines!(ax1, sol.ks .* ρᵢ, real.(ωmax)./ ωₙ)
    scatterlines!(ax2, sol.ks .* ρᵢ, imag.(ωmax)./ ωₙ)
    fig
end
```