#set page(
  paper: "us-letter",
  margin: (x: 0.72in, y: 0.62in),
)

#set par(justify: true, leading: 0.58em)
#set heading(numbering: "1.", outlined: true)
#show heading.where(level: 1): it => block(
  above: 18pt,
  below: 8pt,
  text(size: 15pt, weight: "bold", fill: rgb("#0f172a"), it),
)

#show raw.where(block: true): it => block(
  fill: rgb("#f1f5f9"),
  radius: 3pt,
  inset: 8pt,
  it,
)

#let blue = rgb("#2563eb")
#let cyan = rgb("#0891b2")
#let purple = rgb("#7c3aed")
#let amber = rgb("#d97706")
#let red = rgb("#dc2626")
#let slate = rgb("#64748b")
#let ok = rgb("#166534")

#let chartbar(label, value, maximum, shown, fill: blue) = {
  grid(
    columns: (100pt, 232pt, 72pt),
    gutter: 8pt,
    align: horizon,
    text(size: 9pt, label),
    box(
      width: value / maximum * 226pt,
      height: 11pt,
      fill: fill,
      radius: 2pt,
    ),
    align(right, text(size: 9pt, shown)),
  )
  v(5pt)
}

#let results = toml("results.toml")
#let sharp = toml("sharp_edge.toml")

#align(center)[
  #text(size: 23pt, weight: "bold", fill: rgb("#0f172a"))[
    Arbitrary-Precision BOHH Assembly
  ]
  #v(4pt)
  #text(size: 12pt, fill: slate)[
    Float64, Double64, and BigFloat report
  ]
]

#v(12pt)
#block(fill: rgb("#eff6ff"), stroke: 0.6pt + rgb("#bfdbfe"), radius: 4pt, inset: 11pt)[
  *Finding.* Typed matrix assembly now exposes precision loss in a non-smooth
  Hermite case: relative susceptibility-block error falls from
  #raw(results.errors.Float64) with `Float64` to
  #raw(results.errors.Double64) with `Double64`, and
  #raw(results.errors.BigFloat128) with 128-bit `BigFloat`. On a sharp-edge
  solution comparison, one physical-window root moves by
  #raw(sharp.Float64_error) `Omega_p` in `Float64`, while `Double64` reduces
  that displacement to #raw(sharp.Double64_error) `Omega_p`. Loading
  `GenericLinearAlgebra.jl` enables extended-precision solution.
]

= Question Answered

For Hermite-Hermite (`BOHH`), arithmetic follows the promoted types of
`B0`, `kx`, and `kz`:

```julia
M = dispersion_matrix(species, T(B0), T(kx), T(kz), BOHH;
    N = 2, J = 8)
```

#table(
  columns: (1.38fr, 0.72fr, 2.5fr),
  inset: (x: 6pt, y: 5pt),
  stroke: 0.35pt + rgb("#cbd5e1"),
  fill: (_, y) => if y == 0 { rgb("#e2e8f0") },
  table.header([*Type / package*], [*Assembles?*], [*Status in this worktree*]),
  [`Float64`/`Float32`], [#text(fill: ok)[Yes]], [Default, LAPACK fast path.],
  [`DoubleFloats.Double64`],
  [#text(fill: ok)[Yes]],
  [Finite `Matrix{Complex{Double64}}`; plateau accuracy measured below.],

  [`BigFloat`], [#text(fill: ok)[Yes]], [Finite `Matrix{Complex{BigFloat}}`; tested at 128 and 256 bits.],

  [`MultiFloats.Float64x2`],
  [#text(fill: red)[No]],
  [Blocked: `SpecialFunctions` lacks `erfc`/`besselj` for `Float64x2`. No fallback in PlasmaBO core.],
)

`Double64` support is selected through numeric input types; no required
dependency is added to PlasmaBO. Loading `SpecialFunctions.jl` (and
`DoubleFloats.jl`) activates native special-function dispatch.

== Important boundary

`dispersion_matrix(...)` uses the promoted wave-parameter type. Loading
`GenericLinearAlgebra.jl` enables `solve(...)` for extended scalar types
by extending `LinearAlgebra.eigvals!`; the native `Float64`/`Float32` solve
retains its LAPACK path with no added dependency.

== Selecting precision

`BOHH` is not parameterized by arithmetic type. Supplying extended-precision
wave parameters constructs an extended-precision matrix; loading
`GenericLinearAlgebra.jl` activates the corresponding `solve`.

```julia
using PlasmaBO, DoubleFloats, GenericLinearAlgebra, SpecialFunctions

species = Maxwellian(:e, 1.0, 1.0)
N, J = 1, 4

M64  = dispersion_matrix(species, 1.0, 0.2, 0.1, BOHH; N, J)
Md   = dispersion_matrix(species, Double64(1), Double64(0.2), Double64(0.1), BOHH; N, J)
Mbig = setprecision(128) do
    dispersion_matrix(species, big"1", big"0.2", big"0.1", BOHH; N, J)
end

ω64  = solve(species, 1.0, 0.2, 0.1, BOHH; N, J)         # LAPACK
ωbig = setprecision(128) do
    solve(species, big"1", big"0.2", big"0.1", BOHH; N, J) # GenericLinearAlgebra
end
```

= Regression Case: Non-Smooth Plateau

This test targets Xie's stated failure mode: round-off accumulation becomes
visible for non-smooth velocity distributions and high Hermite order. Input
distribution is a parallel top-hat plateau with smooth perpendicular profile:

```julia
vz = collect(range(-6, 6; length = 1201))
vx = collect(range(0, 6; length = 201))
fv = [abs(z) <= 1 ? exp(-x^2) : 0.0 for z in vz, x in vx]
alm = hermite_expansion(fv, vz, vx, 1.0, 1.0; Nz = 32, Nx = 0).alm

# Assemble BOHH with N = 1, J = 4 and compare its 3 x 3
# accumulated susceptibility block against BigFloat(256).
```

The Hermite projection `alm` is kept identical across types; comparison
therefore isolates loss during eigenproblem assembly, not projection changes.

== Accuracy result

#table(
  columns: (1.45fr, 1.25fr),
  inset: (x: 6pt, y: 5pt),
  stroke: 0.35pt + rgb("#cbd5e1"),
  fill: (_, y) => if y == 0 { rgb("#e2e8f0") },
  table.header([*Assembly arithmetic*], [*Relative block error*]),
  [`Float64`], raw(results.errors.Float64),
  [`Double64`], raw(results.errors.Double64),
  [`BigFloat` (128 bit)], raw(results.errors.BigFloat128),
)

`Double64` supplies about 15 additional reliable decimal digits in this
accumulation. 128-bit `BigFloat` agrees with the 256-bit reference at roughly
37 decimal digits on this measured plateau case.

This is not a general 128-bit accuracy guarantee: generic non-`Float64`
quadrature currently defaults to `sqrt(eps(one(T)))` tolerance. Cases whose
error is dominated by numerical integration may realize fewer digits.

= Sharp-Edge Solution Difference

For solution comparison, an Astfalk-like firehose proton distribution was
truncated with the hard parallel cutoff
`abs(sqrt(2) * v_parallel / vtz) <= 2.35`, then solved at
`Nz = 128`, `N = 2`, and `J = 4`.
`GenericLinearAlgebra.jl` was used on `Double64` and `BigFloat` matrices, while LAPACK
solved the `Float64` matrix. Roots inside
`abs(real(omega / Omega_p)) < 5` and `abs(imag(omega / Omega_p)) < 2` were
matched to a 192-bit `BigFloat` reference by nearest complex distance.
The reference root with largest `Float64` displacement:

#raw(sharp.reference)

#table(
  columns: (1.4fr, 1.25fr),
  inset: (x: 6pt, y: 5pt),
  stroke: 0.35pt + rgb("#cbd5e1"),
  fill: (_, y) => if y == 0 { rgb("#e2e8f0") },
  table.header([*Solved arithmetic*], [*Error / `Omega_p`*]),
  [`Float64` / LAPACK], raw(sharp.Float64_error),
  [`Double64` / generic solve], raw(sharp.Double64_error),
  [`BigFloat` 128 / generic solve], raw(sharp.BigFloat128_error),
)

`GenericLinearAlgebra.jl` was selected after comparison with
`GenericSchur.jl` on the same 132 by 132 sharp-edge `BigFloat(128)` matrix:
median eigensolve was `5.89 s` versus `5.08 s`, but maximum normalized
matched-root error against 192-bit output was `1.98e-19` versus `5.74e-15`.

This is direct solution-level evidence: with the sharp edge, `Float64` fails
to retain this damped branch in its reference location, while extended
arithmetic tracks it. It remains a numerical stress case, not a validation of
the modified distribution as a physical instability model.

#figure(
  image("arbitrary_precision_sharp_edge_distribution.svg", width: 82%),
  caption: [Sharp-edge distribution used for solution comparison: logarithmic density map and integrated parallel lineout.],
)

#figure(
  image("arbitrary_precision_solution_comparison.svg", width: 100%),
  caption: [Sharp-edge solution overlay, zoom around shifted damped branch, and matched-root error by arithmetic type. Frequencies are normalized by proton cyclotron frequency.],
)

= Performance

Measurement uses `Chairmarks.jl` v1.3.1, ten samples with `evals = 1` after
warm-up. Target is the same plateau `dispersion_matrix` (assembly) and full
`solve` (assembly + eigvals). Numbers below are read from
`experiments/arbitrary_precision/results.toml`, regenerated on each
benchmark run.

== Assembly

#table(
  columns: (1.25fr, 1fr, 0.95fr, 1.15fr, 0.9fr),
  inset: (x: 5pt, y: 5pt),
  stroke: 0.35pt + rgb("#cbd5e1"),
  fill: (_, y) => if y == 0 { rgb("#e2e8f0") },
  table.header([*Arithmetic*], [*Median (ms)*], [*Slowdown*], [*Allocations*], [*KiB*]),
  [`Float64`],
  [#str(results.assembly.Float64.median_ms)],
  [#str(results.assembly.Float64.slowdown)x],
  [#str(results.assembly.Float64.allocs)],
  [#str(results.assembly.Float64.kib)],

  [`Double64`],
  [#str(results.assembly.Double64.median_ms)],
  [#str(results.assembly.Double64.slowdown)x],
  [#str(results.assembly.Double64.allocs)],
  [#str(results.assembly.Double64.kib)],

  [`BigFloat` (128)],
  [#str(results.assembly.BigFloat128.median_ms)],
  [#str(results.assembly.BigFloat128.slowdown)x],
  [#str(results.assembly.BigFloat128.allocs)],
  [#str(results.assembly.BigFloat128.kib)],
)

== Solve (assembly + eigvals)

#table(
  columns: (1.25fr, 1fr, 0.95fr, 1.15fr, 0.9fr),
  inset: (x: 5pt, y: 5pt),
  stroke: 0.35pt + rgb("#cbd5e1"),
  fill: (_, y) => if y == 0 { rgb("#e2e8f0") },
  table.header([*Arithmetic*], [*Median (ms)*], [*Slowdown*], [*Allocations*], [*KiB*]),
  [`Float64` / LAPACK],
  [#str(results.solve.Float64.median_ms)],
  [#str(results.solve.Float64.slowdown)x],
  [#str(results.solve.Float64.allocs)],
  [#str(results.solve.Float64.kib)],

  [`Double64` / GLA],
  [#str(results.solve.Double64.median_ms)],
  [#str(results.solve.Double64.slowdown)x],
  [#str(results.solve.Double64.allocs)],
  [#str(results.solve.Double64.kib)],

  [`BigFloat` (128) / GLA],
  [#str(results.solve.BigFloat128.median_ms)],
  [#str(results.solve.BigFloat128.slowdown)x],
  [#str(results.solve.BigFloat128.allocs)],
  [#str(results.solve.BigFloat128.kib)],
)

These costs are PlasmaBO implementation results, not general package
performance claims.

= Limitations

- *J-pole tables*: Xie's published coefficients are finite-precision. Extended
  assembly cannot recover precision absent from those inputs.
- *Generic quadrature tolerance*: `eps(one(T))` per arithmetic. Cases dominated
  by numerical integration may realise fewer digits than the assembly itself.
- *`MultiFloats.Float64x2`*: not supported. `SpecialFunctions` lacks
  `erfc`/`besselj` methods for this type, and PlasmaBO has no generic
  fallback.
- *`BOPBK` eigvals*: assembly works in all supported types, but
  `GenericLinearAlgebra` may hit its Hessenberg iteration limit on the
  ill-conditioned bi-kappa matrix at extended precision. `Float64` LAPACK is
  the validated solve path; `Double64`/`BigFloat` are matrix-construction
  only until the solver is hardened.

= Sources

- #link("https://github.com/JuliaMath/DoubleFloats.jl")[DoubleFloats.jl]:
  `Double64` extended-precision type and package capabilities.
- #link("https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl")[GenericLinearAlgebra.jl]:
  generic dense eigensolution used by the optional `solve` extension.

#pagebreak()
= Appendix A: Sharp-Edge Solution Figure Source

Runnable file: `experiments/arbitrary_precision/precision_solution_plot.jl`.

#{
  set text(size: 7pt)
  raw(read("precision_solution_plot.jl"), lang: "julia", block: true)
}
