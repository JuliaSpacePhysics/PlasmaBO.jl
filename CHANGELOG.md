# Changelog

## [Unreleased]

### Changed

- **Breaking**: plasmas are described with [PlasmaPopulations.jl](https://github.com/JuliaSpacePhysics/PlasmaPopulations.jl) (its main types re-exported). `Maxwellian`, `BiKappa`, `BiKappa2`, `FluidSpecies` and the ChargedParticles/Unitful dependencies are removed; use `Population(particle, BiMaxwellian/GaussianRing/ProductBiKappa/BiKappa(...); n)`. Drifts and ring speeds are in m/s (or Unitful).
- **Breaking**: `HHSolverParam` is internal; give precomputed Hermite coefficients as `HermiteVDF(; w_perp, w_para, alm, vd_para, vr_perp)`. `BOHH` expands any non-Maxwellian model with `HermiteVDF(model)`.
- **Breaking**: kappa models are expanded about their thermal speeds instead of their kappa widths θ, so `BOHH` results for them change; pass `HermiteVDF(model; w_perp, w_para)` for the former reference.
- **Breaking**: `BOHH` expands drifting non-Maxwellian models about their drift; the expansion used to ignore it, so results change.
- **Breaking**: wavevector and algorithm instance are explicit: `solve(species, B0, kx, kz; N)` → `solve(plasma, Wavevector(kx, kz), BOHH(; N))`, `solve(species, B0, ks, θs, alg)` → `solve(plasma, AngleSweep(; k = ks, theta = θs), alg())`. `dispersion_matrix` likewise.
- **Breaking**: fluid closure moved to the solver, `BOFluid(; gamma_perp, gamma_para)`.
- **Breaking**: `hermite_expansion` is public, not exported, and takes vector grids: `hermite_expansion(fv, v_para, v_perp, w_para, w_perp; N_para, N_perp, vd_para, vr_perp, atol)` replaces the keywords `Nz`, `Nx`, `dz`, `dx`. `gen_fv2d` is internal; tabulate through `HermiteVDF(GridVDF(v_perp, v_para, f))`.
- **Breaking**: `BranchPoint`, `SurfaceBranchPoint` are public, not exported; `track(sol, (k, ω))` needs neither.
- **Breaking**: `polarization_ratio(v, k::Wavevector)` and `handedness(v, ω, k::Wavevector; threshold)` replace the `kx`, `kz` arguments and keywords and the defaults to parallel propagation and `ω = 1`.
- **Breaking**: `DispersionSolution` fields are `k`, `theta`, `omega` instead of `ks`, `θs`, `ωs`.
- **Breaking**: `B0` is a magnitude or a vector with `ẑ` along it; a negative scalar no longer reverses the gyration sense.
- **Breaking**: `PlasmaBO.Constants` and the names `q`, `qe`, `me`, `mp`, `kb`, `c0`, `ε0`, `μ0` are removed; use `PlasmaPopulations.E_SI`, `ME_SI`, `MP_SI`, `KB_SI`, `C_SI`, `EPS0_SI`, `MU0_SI`.
- **Breaking**: `plasma_frequency(q, n, m)` and `Debye_length` are removed; use `PlasmaPopulations.plasma_frequency(particle, n)` and `PlasmaPopulations.debye_length(particle, n, T)`.

### Added

- Signed-density components: a `Population` may have negative `n`, enabling exact mixture-of-envelopes representation of subtracted-Maxwellian structures.
- `atol` keyword in `hermite_expansion`: thresholds grid-quadrature noise in high-order coefficients.

## [0.2.0] - 2025-12-26

### Added

- Uniform interface for multiple solver algorithms.
- Add PBK solver.

### Changed

- **Breaking**: Replace `track_dispersion_branch(es)` with `track()`, replace `solve_kinetic_dispersion` with `solve()`. 