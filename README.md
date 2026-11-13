# PlasmaBO

[![DOI](https://zenodo.org/badge/1120953450.svg)](https://doi.org/10.5281/zenodo.18058843)

## Quickstart

PlasmaBO finds every wave mode at once as the eigenvalues of a matrix dispersion relation, no initial guess required. Plasmas are described with [PlasmaPopulations.jl](https://github.com/JuliaPlasma/PlasmaPopulations.jl); problems, solutions and plotting follow [PlasmaDispersionBase.jl](https://github.com/JuliaPlasma/PlasmaDispersionBase.jl#readme).

```julia
using Pkg; Pkg.add("PlasmaBO")
using PlasmaBO

B0 = 100e-9                       # background field [T]
n, T = 8.7e6, 2.857e-3            # density [m^-3], temperature [eV]
p, e = Proton(), Electron()
plasma = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)

# Every mode at one wavevector (k⊥, k∥) [m^-1], in rad/s
k = Wavevector(0.0, 1e-6)
ωs = solve(DispersionProblem(plasma, k), BOHH()).omega

# Every root in an ω box [rad/s] over |k| and propagation angle θ, joined into branches
sweep = AngleSweep(; k = 10 .^ range(-7, -4, length = 50), theta = deg2rad.(0:15:90))
sol = solve(DispersionProblem(plasma, (-1e3 - 1e2im, 1e3 + 1e2im), sweep), BOHH())
filter(isgrowing, link(sol))

# One branch through a seed at the wavevector of the sweep nearest `k`:
# here the whistler wave, among `ωs`
branch = solve(DispersionProblem(plasma, Seed(0.77, k), sweep), BOHH())
branch.omega              # shaped like the sweep's (k, θ) grid
```

`ω = ωr + iγ`: `γ > 0` is an instability, `γ < 0` is damping.

The eigenvalues carry the J-pole and truncation errors of the expansion, and some are spurious.
[VlasovMaxwellDispersion.jl](https://github.com/JuliaSpacePhysics/VlasovMaxwellDispersion.jl) polishes
them on the exact hot-plasma determinant, each root from its own value:

```julia
using VlasovMaxwellDispersion
exact = refine(branch, Muller(); gate = 0.05)   # or a survey; a root moving further than 0.05 Ω_ref is NaN
```

VlasovMaxwellDispersion has no `HermiteVDF`: polish such a plasma's roots on the base models,
`refine(DispersionProblem(model_plasma, sweep), branch.omega, Muller())`.

The default Float64 eigensolver also loses digits, more at larger `J`, and depends on the plasma
(a kappa through `HermiteVDF` loses far more than a Maxwellian at the same `J`).
`refine(sol, BOHH(; N, J); gate)` recovers them: from each root, Newton on `BOHH`'s own rational
relation. A single-wavevector `Seed` and `Continuation(BOHH())` use the same Newton from their guess.

| solver | distributions |
|---|---|
| `BOHH(; N, J)` | using Hermite–Hermite basis; `BiMaxwellian` and `GaussianRing` exactly; other model through `HermiteVDF(model)`; `HermiteVDF(; w_perp, w_para, alm, vd_para, vr_perp)` for given coefficients; harmonics whose coupling is below roundoff are dropped from the matrix, so a generous `N` costs little |
| `BOPBK(; N)` | `ProductBiKappa` with integer `kappa_para` in `1:24` |
| `BOFluid(; gamma_perp, gamma_para)` | multi-fluid solver using population density, drift `vd_para` and `vth_perp_para` |

| name | meaning |
|---|---|
| `Wavevector`, `DispersionProblem`, `solve`, `track`, `Plasma`, `Population`, … | PlasmaDispersionBase's API re-exported |
| `dispersion_matrix(plasma, k::Wavevector, alg; ref)` | the linear operator, eigenvalues in the units of `plasma`; its eigenvectors end in the fields `(E, B)` with `c = 1` |
| `HermiteVDF(model; w_perp, w_para, vd_para, vr_perp, N_para, N_perp, atol)` | expansion of a gyrotropic model, about its thermal speeds and drift by default |
| `polarization_ratio(v, k)` | `Ey/E⊥` of eigenvector `v` |
| `handedness(v, ω, k; threshold)` | `:R`, `:L` or `:linear` |

See [documentation](https://JuliaSpacePhysics.github.io/PlasmaBO.jl/dev/) for worked examples: cold plasma (kinetic vs fluid), ring-beam, firehose, ion-beam, and dispersion-surface tracking.

## Experiments

- [`arbitrary_precision/`](experiments/arbitrary_precision/) — `BOHH`/`BOFluid`/`BOPBK` assembly and solve at `Float64`, `Double64`, and `BigFloat`.
- [`benchmark/botian/`](benchmark/botian/) — speed and root agreement against [BO-T](https://github.com/FusionAlpha/botian) (Python).

## Elsewhere

- [hsxie/BO-Arbitrary](https://github.com/hsxie/boarbitrary/tree/main): Extension of the kinetic electromagnetic magnetized dispersion relation solver [PDRK](https://github.com/hsxie/pdrk)/[BO](https://github.com/hsxie/bo) to arbitrary distributions (MATLAB)
  - [liangwang0734/xenon](https://github.com/liangwang0734/xenon): A matriX-based dispErsioN relatiON solver (Python). Very limited in functionality (only electrostatic for kinetic BiMaxwellian plasmas).
- [LinearMaxwellVlasov.jl](https://github.com/jwscook/LinearMaxwellVlasov.jl): A Julia implementation of linear Maxwell-Vlasov solvers
- [danielver02/ALPS](https://github.com/danielver02/ALPS): The Arbitrary Linear Plasma Solver that solves the Vlasov-Maxwell dispersion relation in hot (even relativistic) magnetised plasma (Fortran)
- [pastfalk/LEOPARD](https://github.com/pastfalk/LEOPARD): Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions (Fortran)
- [Drakicy/MPDES](https://github.com/Drakicy/MPDES): Magnetized Plasma Dispersion Equation Solver (MPDES) (MATLAB)
