"""
    RelativisticLongitudinal

Experimental: dispersion of the longitudinal (Langmuir) mode in an unmagnetized,
isotropic Maxwell–Jüttner plasma.

Provides three independent solvers for the same problem so they can be cross-checked:

  * [`solve_relativistic_direct`](@ref) — direct numerical integration of the
    relativistic dielectric kernel; treated as ground truth.
  * [`solve_relativistic_matrix`](@ref) — BO-style matrix eigenvalue solver
    via Padé approximation of Lerche's F-function and Gauss–Laguerre quadrature
    in γ. Returns all roots at once, no initial guess. Valid for super-luminal
    modes (no Landau damping captured by the real-pole Padé).
  * [`solve_maxwellian_matrix`](@ref) — standard non-relativistic matrix
    method (J-pole expansion of the Z-function) for the same problem.
    The fully classical baseline.

Conventions: ω̃ ≡ ω/ω_p, k̃ ≡ k λ_D, μ ≡ mc²/T (T = electron temperature).
The classical limit is μ → ∞ (T ≪ mc²).
"""
module RelativisticLongitudinal

using LinearAlgebra
using QuadGK: quadgk
using Bessels: besselkx
using ..PlasmaBO: get_jpole_coefficients

export solve_relativistic_direct, solve_relativistic_matrix, solve_maxwellian_matrix
export pade_F, gauss_laguerre, juttner_dielectric

# ===========================================================================
# Padé approximation of Lerche's F-function
# ===========================================================================

"""
    pade_F(J::Int) -> (a::Vector, r::Vector)

[J-1/J] Padé approximant of F(s) = 1 − (s/2)·log((s+1)/(s-1)) in the variable
t = 1/s², matching the first 2J terms of the asymptotic expansion
F(t) = -Σ_{k≥1} t^k/(2k+1).

Returns coefficients (a, r) such that

    F(s) ≈ Σⱼ aⱼ / (s² − rⱼ),   rⱼ = sⱼ²

For all J ≥ 1, the poles rⱼ are real and lie in (0, 1) — they approximate
the s ∈ [−1, 1] log branch cut. Consequently the approximation is accurate
only for |s| ≳ 1 (super-luminal phase velocities).

Padé saturates around J = 4–6; higher J yields no improvement until quadruple
precision is used in the Hankel solve.
"""
function pade_F(J::Int)
    Fcoef(k) = k == 0 ? 0.0 : -1 / (2k + 1)
    A = [Fcoef(J + r - j) for r in 1:J, j in 1:J]
    b = [-Fcoef(J + r) for r in 1:J]
    c = A \ b
    # Denominator Q(t) = 1 + c₁t + … + c_J t^J → roots via companion matrix
    poly = [1.0; c] ./ c[end]
    n = J
    C = zeros(n, n)
    for i in 1:n-1
        C[i+1, i] = 1.0
    end
    C[:, n] = -poly[1:n]
    tj = eigvals(C)
    rs = 1 ./ tj
    # Residues via Vandermonde from first J moments
    V = [rs[j]^(k - 1) for k in 1:J, j in 1:J]
    a = V \ [-1 / (2k + 1) for k in 1:J]
    return a, rs
end

# ===========================================================================
# Gauss–Laguerre quadrature (Golub–Welsch)
# ===========================================================================

"""
    gauss_laguerre(n::Int) -> (x, w)

Gauss–Laguerre nodes and weights for ∫₀^∞ f(x) e⁻ˣ dx ≈ Σ wᵢ f(xᵢ),
computed by eigen-decomposition of the Laguerre Jacobi tridiagonal.
"""
function gauss_laguerre(n::Int)
    α = Float64[2k + 1 for k in 0:n-1]
    β = Float64[k for k in 1:n-1]
    F = eigen(SymTridiagonal(α, β))
    return F.values, F.vectors[1, :] .^ 2
end

# ===========================================================================
# Direct dielectric for unmagnetized Maxwell–Jüttner background
# ===========================================================================

"""
    juttner_dielectric(ω̃, k̃, μ; landau_branch=true) -> ComplexF64

Longitudinal dielectric ε_l(ω, k) for a single-species electron plasma with
isotropic Maxwell–Jüttner equilibrium of inverse temperature μ = mc²/T.

When `landau_branch=true` (default), the Landau prescription ω → ω + i 0⁺ is
applied via principal-value split + closed-form residue (giving Im ε > 0 in
the resonance region β > u). When `false`, only the principal-value real part
is computed.

This is the **direct integration ground truth** used by
[`solve_relativistic_direct`](@ref).
"""
function juttner_dielectric(ω̃::Real, k̃::Real, μ::Real; landau_branch::Bool = true)
    u = ω̃ / (k̃ * sqrt(μ))
    K2x = besselkx(2, μ)
    γ_res = u < 1 ? 1 / sqrt(1 - u^2) : Inf
    # Real part — principal value
    f(γ) = begin
        β = sqrt(γ^2 - 1) / γ
        (γ^2 - 1) * exp(-(γ - 1) * μ) * (2 / β - (u / β^2) * log(abs((u + β) / (u - β))))
    end
    if isinf(γ_res)
        J_real, _ = quadgk(f, 1.0, Inf; rtol = 1e-10)
    else
        δ = min(0.01, 0.05 * (γ_res - 1))
        J1, _ = quadgk(f, 1.0, γ_res - δ; rtol = 1e-10)
        J2, _ = quadgk(f, γ_res + δ, Inf; rtol = 1e-10)
        J_real = J1 + J2
    end
    re = 1 + (μ / (2 * k̃^2 * K2x)) * J_real

    if !landau_branch || isinf(γ_res)
        return complex(re, 0.0)
    end

    g(γ) = begin
        β = sqrt(γ^2 - 1) / γ
        (γ^2 - 1) * exp(-(γ - 1) * μ) * (u / β^2)
    end
    Ir, _ = quadgk(g, γ_res, Inf; rtol = 1e-10)
    im_part = +(μ / (2 * k̃^2 * K2x)) * π * Ir
    return complex(re, im_part)
end

# ===========================================================================
# Method 1: direct integration
# ===========================================================================

"""
    omega0_asymptotic(μ) -> Real

Relativistic Langmuir frequency at k → 0, ω̃₀² = μ · ⟨β²⟩_J / 3, where
⟨β²⟩_J is the Jüttner expectation of v²/c². Used as seed for root-finding.

Limits: ω̃₀ → 1 as μ → ∞ (classical), ω̃₀ → √(μ/3) → 0 as μ → 0 (ultra-rel).
"""
function omega0_asymptotic(μ::Real)
    val, _ = quadgk(γ -> ((γ^2 - 1)^1.5 / γ) * exp(-(γ - 1) * μ), 1.0, Inf; rtol = 1e-12)
    β² = (μ / besselkx(2, μ)) * val
    return sqrt(μ * β² / 3)
end

"""
    solve_relativistic_direct(k̃, μ; ωlow=nothing, ωhi=nothing, n=500) -> (ωr, γ)

Find the Langmuir mode by bisecting Re ε_l(ω) = 0, then estimate the damping
rate γ from the weak-damping formula γ ≈ −Im ε_l / ∂_ω Re ε_l.

By default scans an ω̃ window centered on the asymptotic k → 0 frequency
ω̃₀(μ) = √(μ ⟨β²⟩_J / 3) (see [`omega0_asymptotic`](@ref)). This skips
spurious roots that the principal-value integrand can show at small ω̃.

Ground-truth solver for cross-validation. Returns one root per call.
"""
function solve_relativistic_direct(
    k̃::Real, μ::Real;
    ωlow::Union{Real,Nothing} = nothing,
    ωhi::Union{Real,Nothing} = nothing,
    n::Int = 500,
)
    ω0 = omega0_asymptotic(μ)
    # Window must span from below ω₀ to above the Bohm-Gross-shifted root.
    # Take a generous upper bound: ω₀ + 3·k̃ (approx leading thermal shift).
    ωlow_eff = something(ωlow, max(0.5 * ω0, 0.05))
    ωhi_eff = something(ωhi, max(ω0 + 3.0 * k̃ + 0.2, 1.5))
    fr(ω) = real(juttner_dielectric(ω, k̃, μ))
    # Scan top-down: the *highest-frequency* sign-change is the Langmuir branch.
    # (A second sub-Langmuir root can appear at large k̃ from the principal-value
    # integral; it is a separate physical mode, not the wave we want here.)
    ωs = range(ωhi_eff, ωlow_eff, length = n)
    last = fr(ωs[1])
    a, b = NaN, NaN
    for i in 2:n
        cur = fr(ωs[i])
        if isfinite(cur * last) && cur * last < 0
            a, b = ωs[i], ωs[i-1]   # restore lo < hi ordering
            break
        end
        last = cur
    end
    isnan(a) && return (NaN, NaN)
    for _ in 1:80
        m = 0.5 * (a + b)
        fr(a) * fr(m) < 0 ? (b = m) : (a = m)
        b - a < 1e-10 && break
    end
    ωr = 0.5 * (a + b)
    h = 1e-4
    dε = (fr(ωr + h) - fr(ωr - h)) / (2h)
    γ = -imag(juttner_dielectric(ωr, k̃, μ)) / dε
    return (ωr, γ)
end

export omega0_asymptotic

# ===========================================================================
# Method 2: relativistic matrix solver via Padé-of-F
# ===========================================================================

"""
    solve_relativistic_matrix(k̃, μ; J=6, N=20) -> Vector{ComplexF64}

Returns *all* eigenvalues ω̃ of the relativistic dispersion relation built by:

  1. [J-1/J] Padé of Lerche's F(s) → J poles in s² (real, in (0,1)).
  2. Gauss–Laguerre quadrature with N nodes in γ ∈ [1, ∞).
  3. Companion-matrix linearisation (BO trick).

Matrix dimension is 2·J·N + 1 (≈ 241 for default J=6, N=20). Per-call cost
is dominated by LAPACK `eigvals` on this dense matrix (≈ 40 ms).

The physical Langmuir root is the largest-magnitude real eigenvalue;
others are quadrature-artifact "ballistic" modes clustered near the pole
positions ±β(γ_q)·s_j·k̃√μ.

Damping is NOT captured because all Padé poles are real (see Note in module
docstring). Extension via complex-pole rational approximation (AAA / Carrier–
Krook–Pearson) would recover Landau damping.
"""
function solve_relativistic_matrix(k̃::Real, μ::Real; J::Int = 6, N::Int = 20, with_eigvecs::Bool = false)
    a, rs = pade_F(J)
    xq, wq = gauss_laguerre(N)
    γs = 1 .+ xq ./ μ
    βs = sqrt.(1 .- 1 ./ γs .^ 2)
    K2x = besselkx(2, μ)
    poles = ComplexF64[]
    coefs = ComplexF64[]
    for j in 1:J
        sj = sqrt(complex(rs[j]))
        for q in 1:N
            wq_eff = wq[q] / μ
            φ = βs[q] * sj * k̃ * sqrt(μ)
            Ajq = a[j] * wq_eff * (γs[q]^2 - 1) * μ^2 /
                  (2 * sj * k̃ * sqrt(μ) * K2x)
            push!(poles, φ)
            push!(coefs, Ajq)
            push!(poles, -φ)
            push!(coefs, -Ajq)
        end
    end
    M = length(poles)
    Mat = zeros(ComplexF64, M + 1, M + 1)
    Mat[1, 1] = -sum(coefs)
    @inbounds for p in 1:M
        Mat[1, p+1] = -poles[p]
        Mat[p+1, 1] = coefs[p]
        Mat[p+1, p+1] = poles[p]
    end
    if with_eigvecs
        F = eigen(Mat)
        return F.values, F.vectors
    else
        return eigvals(Mat)
    end
end

# ===========================================================================
# Method 3: non-relativistic Maxwellian matrix (J-pole Z) — baseline
# ===========================================================================

"""
    solve_maxwellian_matrix(k̃; J=8) -> Vector{ComplexF64}

Returns all eigenvalues of the *classical* (non-relativistic) Langmuir
dispersion relation 1 + (1/k̃²)[1 + ζZ(ζ)] = 0, with ζ = ω̃/(√2 k̃),
using the standard J-pole expansion Z(ζ) ≈ Σⱼ bⱼ/(ζ − cⱼ) from
[`PlasmaBO.get_jpole_coefficients`](@ref).

Uses the identity Σbⱼ = -1 to write the dispersion as
ε_l − 1 = (√2/k̃) Σⱼ bⱼ cⱼ /(ω̃ − √2 k̃ cⱼ), giving a (J+1)×(J+1) companion
matrix — tiny and instantaneous. This is the natural baseline against
which the relativistic solvers are compared at large μ.

Note: this solver has no μ dependence — that is exactly the point of the
comparison. The classical limit assumes all particles move slower than
the speed of light and the dispersion depends only on k̃.
"""
function solve_maxwellian_matrix(k̃::Real; J::Int = 8)
    coeffs = get_jpole_coefficients(J)
    bj = coeffs.bzj
    cj = coeffs.czj
    M = length(bj)
    # ε_l − 1 = Σⱼ Bⱼ/(ω̃ − pⱼ),  Bⱼ = (√2/k̃) bⱼ cⱼ,  pⱼ = √2 k̃ cⱼ
    # ε_l = 1 + (1/k̃²)[1 + ζZ] with ζ = ω̃/(√2 k̃) and Σbⱼ = -1 gives
    #   ε_l − 1 = (1/k̃²) Σⱼ bⱼ cⱼ /(ζ − cⱼ) = (√2/k̃) Σⱼ bⱼ cⱼ /(ω̃ − √2 k̃ cⱼ)
    poles = ComplexF64[sqrt(2) * k̃ * cj[j] for j in 1:M]
    coefs = ComplexF64[(sqrt(2) / k̃) * bj[j] * cj[j] for j in 1:M]
    Mat = zeros(ComplexF64, M + 1, M + 1)
    Mat[1, 1] = -sum(coefs)
    @inbounds for p in 1:M
        Mat[1, p+1] = -poles[p]
        Mat[p+1, 1] = coefs[p]
        Mat[p+1, p+1] = poles[p]
    end
    return eigvals(Mat)
end

# ===========================================================================
# Helper: pick the physical Langmuir branch from a vector of eigenvalues
# ===========================================================================

"""
    langmuir_root(eigs; near=nothing, ωmin=0.05) -> ComplexF64 | nothing

Pick a Langmuir-like root from a vector of dispersion eigenvalues.

  * If `near` is given: returns the eigenvalue closest to it (in real part)
    among those with positive real part and `|imag| < 0.1`.
  * If `near` is `nothing`: returns the largest-real-part eigenvalue
    above `ωmin` — adequate for the non-relativistic matrix method where
    no ballistic eigenvalues are present.

For the relativistic matrix solver, **always pass `near`**: the spectrum
contains O(JN) ballistic eigenvalues clustered along the real axis and the
physical Langmuir mode is one specific eigenvalue interleaved among them.
The natural seed is [`omega0_asymptotic`](@ref) (the k → 0 relativistic
plasma frequency) plus a leading-order thermal shift if available.
"""
function langmuir_root(eigs::AbstractVector{<:Complex};
                       near::Union{Real,Complex,Nothing} = nothing,
                       ωmin::Real = 0.05)
    cand = filter(e -> real(e) > ωmin && abs(imag(e)) < 0.1, eigs)
    isempty(cand) && return nothing
    if near === nothing
        return cand[argmax(real.(cand))]
    else
        return cand[argmin(abs.(real.(cand) .- real(near)))]
    end
end

export langmuir_root

end # module
