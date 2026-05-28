# Arbitrary precision

`BOHH`/`BOFluid`/`BOPBK` matrix assembly at `Float64`, `Double64`, and `BigFloat`.
Solve via LAPACK (`Float32`/`Float64`) or `GenericLinearAlgebra` (extended types).

📄 [Rendered report (PDF, v0.2.2)](https://github.com/JuliaSpacePhysics/PlasmaBO.jl/releases/download/v0.2.2/arbitrary_precision_report.pdf)

## Reproduce

```sh
julia --project=. -e 'include("./arbitrary_precision_report.jl")'
julia --project=. -e 'include("./precision_solution_plot.jl")'
typst compile arbitrary_precision_report.typ
```

The two `.jl` drivers emit `results.toml` and `sharp_edge.toml` plus two plots.
The Typst source consumes them.

## Files

| File | Role |
|---|---|
| `arbitrary_precision_report.typ` | Typst source  |
| `arbitrary_precision_report.jl` | Plateau accuracy + assembly/solve benchmarks → `results.toml` |
| `precision_solution_plot.jl` | Sharp-edge solution comparison + figure → `sharp_edge.toml` + SVGs |
| `besselj_compare.jl` | Bessel back-end comparison probe |
