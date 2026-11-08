# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)

## Quickstart

`solve` returns all wave modes at once as eigenvalues of a matrix dispersion relation — no initial guess required. Describe the plasma with [PlasmaPopulations.jl](https://github.com/JuliaSpacePhysics/PlasmaPopulations.jl) (its main types re-exported), pick a solver (`BOHH` kinetic, Hermite–Hermite basis for arbitrary distributions), and call:

```julia
using Pkg; Pkg.add("PlasmaBO")
using PlasmaBO

B0 = 100e-9                       # background field [T]
n, T = 8.7e6, 2.857e-3            # density [m^-3], temperature [eV]
p, e = Proton(), Electron()
plasma = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)

# One wavevector (k⊥, k∥) [m^-1]: Vector{ComplexF64} of all ω roots [rad/s]
k = Wavevector(0.0, 1e-6)
ωs = solve(plasma, k, BOHH())
ωs_fluid = solve(plasma, k, BOFluid())

# Scan over |k| and propagation angle θ
k = 10 .^ range(-7, -4, length = 50)
theta = deg2rad.(0:15:90)
sol = solve(plasma, AngleSweep(; k, theta), BOHH())   # sol.omega[ik, itheta] :: Vector{ComplexF64}
```

`ω = ωr + iγ`: `γ > 0` is an instability, `γ < 0` is damping.

| solver | distributions |
|---|---|
| `BOHH(; N, J)` | using Hermite–Hermite basis; `BiMaxwellian` and `GaussianRing` exactly; other model through `HermiteVDF(model)`; `HermiteVDF(; w_perp, w_para, alm, vd_para, vr_perp)` for given coefficients |
| `BOPBK(; N)` | `ProductBiKappa` with integer `kappa_para` |
| `BOFluid(; gamma_perp, gamma_para)` | multi-fluid solver using population density, drift `vd_para` and `vth_perp_para` |

| name | meaning |
|---|---|
| `Particle`, `Electron`, `Positron`, `Proton`, `Population`, `Plasma`, `BiMaxwellian`, `GaussianRing`, `BiKappa`, `ProductBiKappa`, `Cold`, `GridVDF`, `Wavevector`, `AngleSweep` | re-exported from PlasmaPopulations; `using PlasmaPopulations` for the rest of it |
| `solve(plasma, k::Wavevector, alg)` | `Vector` of all `ω` [rad/s] |
| `solve(plasma, sweep::AngleSweep, alg)` | `DispersionSolution` with fields `k`, `theta`, `omega`, `omega[ik, itheta]` being a `Vector` |
| `dispersion_matrix(plasma, k, alg)` | the linear operator; its eigenvectors end in `(E, B)` |
| `HermiteVDF(model; w_perp, w_para, vd_para, vr_perp, N_para, N_perp, atol)` | expansion of a gyrotropic model, about its thermal speeds and drift by default |
| `track(sol, (k, ω))`, `track(sol, (k, θ, ω))` | follow one branch from a seed: `(ks, ωs)`, or surface `(ks, θs, ωs)` |
| `plot_branches(branches, kn, ωn)`, `plot(sol, kn, ωn)` | with Makie loaded |
| `polarization_ratio(v, k)` | `Ey/E⊥` of eigenvector `v` |
| `handedness(v, ω, k; threshold)` | `:R`, `:L` or `:linear` |
| below | public, not exported |
| `PlasmaBO.hermite_expansion(fv, v_para, v_perp, w_para, w_perp; N_para, N_perp, vd_para, vr_perp, atol)` | coefficients `(; alm, a0lm)` of a table `fv[i_para, i_perp]` on uniform grids |
| `PlasmaBO.BranchPoint(k, ω, locator)`, `PlasmaBO.SurfaceBranchPoint(k, θ, ω, locator)` | seeds for `track` with a custom eigenvalue selector |
| `PlasmaBO.DispersionSolution` | return type of a sweep |

See [documentation](https://JuliaSpacePhysics.github.io/PlasmaBO.jl/dev/) for worked examples: cold plasma (kinetic vs fluid), ring-beam, firehose, ion-beam, and dispersion-surface tracking.

## Experiments

- [`arbitrary_precision/`](experiments/arbitrary_precision/) — `BOHH`/`BOFluid`/`BOPBK` assembly and solve at `Float64`, `Double64`, and `BigFloat`.

## Elsewhere

- [hsxie/BO-Arbitrary](https://github.com/hsxie/boarbitrary/tree/main): Extension of the kinetic electromagnetic magnetized dispersion relation solver [PDRK](https://github.com/hsxie/pdrk)/[BO](https://github.com/hsxie/bo) to arbitrary distributions (MATLAB)
  - [liangwang0734/xenon](https://github.com/liangwang0734/xenon): A matriX-based dispErsioN relatiON solver (Python). Very limited in functionality (only electrostatic for kinetic BiMaxwellian plasmas).
- [LinearMaxwellVlasov.jl](https://github.com/jwscook/LinearMaxwellVlasov.jl): A Julia implementation of linear Maxwell-Vlasov solvers
- [danielver02/ALPS](https://github.com/danielver02/ALPS): The Arbitrary Linear Plasma Solver that solves the Vlasov-Maxwell dispersion relation in hot (even relativistic) magnetised plasma (Fortran)
- [pastfalk/LEOPARD](https://github.com/pastfalk/LEOPARD): Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions (Fortran)
- [Drakicy/MPDES](https://github.com/Drakicy/MPDES): Magnetized Plasma Dispersion Equation Solver (MPDES) (MATLAB)
