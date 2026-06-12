# Changelog

## [Unreleased]

### Added

- Signed-density components in `BOHH`: `Maxwellian` (and the direct `HHSolverParam` constructor) accept negative `n`, enabling exact mixture-of-envelopes representation of subtracted-Maxwellian structures.
- `atol` keyword in `hermite_expansion`: thresholds grid-quadrature noise in high-order coefficients.

## [0.2.0] - 2025-12-26

### Added

- Uniform interface for multiple solver algorithms.
- Add PBK solver.

### Changed

- **Breaking**: Replace `track_dispersion_branch(es)` with `track()`, replace `solve_kinetic_dispersion` with `solve()`. 