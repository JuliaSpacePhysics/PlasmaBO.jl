# Relativistic Langmuir wave — three solvers compared

This page compares three independent solvers for the same problem — the
longitudinal (Langmuir) mode in an unmagnetized, isotropic plasma —
implemented in [`PlasmaBO.RelativisticLongitudinal`](@ref).

| Method | Function | Approach |
|:---|:---|:---|
| **Direct integration** | [`solve_relativistic_direct`](@ref) | `QuadGK` on the relativistic dielectric kernel + bisection. Single root per call, used as ground truth. |
| **Relativistic matrix** | [`solve_relativistic_matrix`](@ref) | BO-style: Padé approximation of Lerche's F-function + Gauss–Laguerre quadrature in γ + companion-matrix eigenvalues. Returns all roots simultaneously. |
| **Non-relativistic matrix** | [`solve_maxwellian_matrix`](@ref) | Classical Maxwellian Langmuir wave via standard J-pole expansion of Z(ζ). Independent of μ. Baseline. |

The three methods should converge in the classical limit μ ≡ mc²/T → ∞.
They will *differ* exactly where relativistic kinematics matter.

## Background — the dielectric

For an isotropic Maxwell–Jüttner background of inverse temperature μ:

```math
\varepsilon_l(\tilde\omega, \tilde k) = 1 + \frac{\mu}{2 \tilde k^2 \, K_2(\mu) e^{\mu}}
  \int_1^{\infty} (\gamma^2-1) \, e^{-(\gamma-1)\mu}
  \left[\frac{2}{\beta} - \frac{u}{\beta^2} \ln\!\frac{u+\beta}{u-\beta}\right] d\gamma
```

with ``u = \tilde\omega / (\tilde k \sqrt\mu)``, ``\beta = \sqrt{\gamma^2-1}/\gamma``,
``\tilde\omega = \omega/\omega_p``, ``\tilde k = k \lambda_D``. The bracketed
factor is ``(2/\beta)\,F(u/\beta)`` where ``F(s) = 1 - (s/2)\ln((s+1)/(s-1))``
is Lerche's F-function.

The k → 0 asymptote is ``\tilde\omega_0^2(\mu) = \mu\,\langle\beta^2\rangle_J/3``,
implemented as [`omega0_asymptotic`](@ref). Limits:
- μ → ∞ (classical): ``\tilde\omega_0 \to 1`` (standard ``\omega_p``).
- μ → 0 (ultra-rel): ``\tilde\omega_0 \to \sqrt{\mu/3} \to 0``.

## How each method works

**Direct integration.** Adaptive `QuadGK` on the principal value of the
γ-integral, with a Landau-residue closed form for the imaginary part. Then
bisect Re ε_l = 0. Bracket is auto-chosen around ``\tilde\omega_0(\mu)`` and
scanned top-down so the highest-frequency root (the Langmuir branch) is
returned even when a sub-Langmuir companion exists.

**Relativistic matrix.** [J−1/J] Padé of F(s) in t = 1/s² gives J real
poles ``s_j^2 \in (0, 1)`` — they approximate the s ∈ [-1, 1] log branch
cut. Substituting into the kernel and using partial fractions:

```math
\varepsilon_l(\tilde\omega) - 1 = \sum_{j,q}
  \left[ \frac{A_{jq}}{\tilde\omega - \varphi_{jq}} - \frac{A_{jq}}{\tilde\omega + \varphi_{jq}} \right],
\quad \varphi_{jq} = \beta(\gamma_q) \, s_j \, \tilde k \sqrt\mu
```

with Gauss–Laguerre nodes ``\{\gamma_q, w_q\}`` on γ ∈ [1, ∞) (e⁻ˣ weighted).
The 2JN simple poles are linearised into a companion matrix of size 2JN+1;
its eigenvalues are *all* roots of ε_l = 0. Default J=6, N=20 → 241×241
dense eigenproblem (≈ 40 ms via LAPACK).

**Non-relativistic matrix.** Standard `Z(ζ) ≈ Σⱼ bⱼ/(ζ − cⱼ)` from
[`PlasmaBO.get_jpole_coefficients`](@ref). With ``\zeta = \tilde\omega/(\sqrt 2 \tilde k)``
and the identity ``\sum b_j = -1``:

```math
\varepsilon_l(\tilde\omega) - 1 = \frac{\sqrt 2}{\tilde k}
  \sum_j \frac{b_j c_j}{\tilde\omega - \sqrt 2 \tilde k c_j}.
```

J=8 J-pole gives a 9×9 matrix — tiny and instantaneous. The classical
Z-function admits complex poles, so this matrix naturally captures Landau
damping (the relativistic Padé-of-F variant does not yet — see Note below).

## Side-by-side comparison

```@example relcomp
using PlasmaBO
using Printf

println("μ      k̃      ω̃₀       direct     rel-matrix    Maxw-matrix")
println(repeat("-", 64))
for μ in [1000.0, 100.0, 10.0, 1.0, 0.3]
    ω0 = omega0_asymptotic(μ)
    for k̃ in [0.05, 0.1, 0.2, 0.3]
        ωdir, _ = solve_relativistic_direct(k̃, μ)
        ωrel = langmuir_root(
            solve_relativistic_matrix(k̃, μ; J=6, N=20);
            near = ω0 + 1.5*k̃^2/ω0,
        )
        ωmax = langmuir_root(solve_maxwellian_matrix(k̃; J=8))
        rrel = ωrel === nothing ? NaN : real(ωrel)
        rmax = ωmax === nothing ? NaN : real(ωmax)
        @printf("%5.1f  %.2f   %.4f   %.5f    %.5f      %.5f\n",
                μ, k̃, ω0, ωdir, rrel, rmax)
    end
end
```

Each method tracks the Langmuir mode differently as μ varies:

  * **At large μ** (classical, T ≪ mc²) all three methods agree.
  * **At moderate μ** (T ~ mc²/10) the two relativistic methods agree
    with each other but disagree with the Maxwellian baseline by O(1) —
    *this gap is exactly the relativistic correction the method is
    designed to capture*.
  * **At very small μ** (ultra-relativistic, T ≳ mc²) the relativistic
    methods give ``\tilde\omega \ll 1`` while the Maxwellian method still
    predicts ``\tilde\omega \approx 1`` — the latter is just wrong.

## Plots

```@example relcomp
using CairoMakie

ks = 0.02:0.01:0.45
μs = [1000.0, 100.0, 10.0, 1.0]

fig = Figure(size = (820, 520), fontsize = 14)
ax = Axis(fig[1,1],
    xlabel = "k λ_D", ylabel = "ω_r / ω_p",
    title = "Langmuir dispersion — three solvers compared")

colors = Makie.wong_colors()

_real_or_nan(x) = x === nothing ? NaN : real(x)

# Maxwellian baseline (μ-independent): single curve
ωmax_curve = Float64[
    _real_or_nan(langmuir_root(solve_maxwellian_matrix(k̃; J=8)))
    for k̃ in ks
]
lines!(ax, ks, ωmax_curve;
       color = :black, linewidth = 2.5, linestyle = :dot,
       label = "Maxwell matrix (no μ dep.)")

# Direct integration + relativistic matrix per μ
for (i, μ) in enumerate(μs)
    ω0 = omega0_asymptotic(μ)
    ωdir = Float64[]
    ωrel = Float64[]
    for k̃ in ks
        wd, _ = solve_relativistic_direct(k̃, μ)
        push!(ωdir, wd)
        push!(ωrel, _real_or_nan(langmuir_root(
            solve_relativistic_matrix(k̃, μ; J=6, N=20);
            near = ω0 + 1.5*k̃^2/ω0,
        )))
    end
    lines!(ax, ks, ωdir; color = colors[i], linewidth = 2.2,
           label = "direct  μ=$(μ)")
    scatter!(ax, ks, ωrel; color = colors[i], marker = :cross, markersize = 9,
             label = "rel-matrix  μ=$(μ)")
end

xlims!(ax, 0, 0.45); ylims!(ax, 0, 1.4)
axislegend(ax; position = :lt, nbanks = 2, framevisible = true)
save("rel_comparison.png", fig; px_per_unit = 2)
fig
```

```@example relcomp
# Error plot: relativistic matrix vs direct
fig2 = Figure(size = (820, 380), fontsize = 14)
ax2 = Axis(fig2[1,1],
    xlabel = "k λ_D", ylabel = "|ω̃_matrix − ω̃_direct|",
    title = "Relativistic matrix vs direct integration (J=6, N=20)",
    yscale = log10)

for (i, μ) in enumerate(μs)
    ω0 = omega0_asymptotic(μ)
    xs = Float64[]; ys = Float64[]
    for k̃ in ks
        wd, _ = solve_relativistic_direct(k̃, μ)
        isnan(wd) && continue
        e = langmuir_root(
            solve_relativistic_matrix(k̃, μ; J=6, N=20);
            near = ω0 + 1.5*k̃^2/ω0,
        )
        e === nothing && continue
        push!(xs, k̃); push!(ys, max(1e-10, abs(real(e) - wd)))
    end
    lines!(ax2, xs, ys; color = colors[i], linewidth = 2,
           label = "μ = $(μ)")
    scatter!(ax2, xs, ys; color = colors[i], markersize = 6)
end
xlims!(ax2, 0, 0.45); ylims!(ax2, 1e-7, 1e-1)
axislegend(ax2; position = :rb, framevisible = true)
save("rel_comparison_err.png", fig2; px_per_unit = 2)
fig2
```

## Reading the comparison

The matrix-vs-direct error stays at ~10⁻⁵ over the range where Landau
damping is absent (u = ω/(kc) > 1, super-luminal). It grows by 2–3 decades
once the wave becomes sub-luminal and Landau damping turns on, because the
real-pole Padé of F cannot represent the log-cut residue.

The Maxwellian matrix dispersion (dotted) ignores μ entirely. The gap
between it and the direct curve is the *relativistic correction*. For μ = 1
(T = mc²) the gap is comparable to ω_p itself — the classical solver is
roughly 100 % off and the relativistic one is necessary.

## When to use which solver

  * **Use the direct solver** when you need γ (Landau damping) in the
    weak-damping regime, or for cross-checking a matrix-method root.
  * **Use the relativistic matrix solver** for ω(k) sweeps and dispersion
    surfaces — it returns all roots at once, no continuation needed.
    Current limitation: real-pole Padé does not capture Landau damping;
    output Im(ω) is zero to numerical precision. Damping requires either
    a complex-pole rational approximation (AAA / Carrier–Krook–Pearson)
    or explicit residue injection.
  * **Use the Maxwellian matrix solver** when μ ≫ 1000 or when you only
    care about the classical limit. It's a 9×9 eigenproblem and is
    essentially free.

## Performance

| Solver | Cost per (k̃, μ) | Notes |
|:---|:---|:---|
| Direct integration | 10–50 ms | bracket scan + bisect; single root |
| Relativistic matrix (J=6, N=20) | ≈ 40 ms | 241×241 dense `eigvals`; all roots |
| Maxwellian matrix (J=8) | ≈ 0.1 ms | 9×9 eigenproblem; all roots |

## API reference

```@docs
PlasmaBO.RelativisticLongitudinal
solve_relativistic_direct
solve_relativistic_matrix
solve_maxwellian_matrix
juttner_dielectric
omega0_asymptotic
langmuir_root
PlasmaBO.RelativisticLongitudinal.pade_F
PlasmaBO.RelativisticLongitudinal.gauss_laguerre
```
