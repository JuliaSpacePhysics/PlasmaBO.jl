# Fill `Γ` for harmonic `n ≥ 0` and all m in 0:(m_max+2) with a single adaptive quadrature over
# y ∈ [max(0, d - ym), d + ym]. With V = (y·nJₙ(ay)/(ay), y·Jₙ′(ay), Jₙ(ay)) and its six products
# Vᵢ·Vⱼ in the order (11, 12, 13, 22, 23, 33), the q-th product gives
#   Γ[2q - 1 + p, m + 1] = ∫ (Vᵢ·Vⱼ) · exp(-(y-d)²) · (y-d)^m · y^p dy,   p ∈ (0, 1).
# nJₙ(ay)/(ay) carries the susceptibility's n·ωc/(k⊥v⊥) into the integrand, regular at k⊥ = 0.
function _compute_integral_matrices!(Γₙ, n, a, d, m_max; rtol = nothing, atol = nothing)
    T = eltype(Γₙ)
    tol = T === Float64 ? T(1.0e-10) : sqrt(eps(one(T)))
    rtol = something(rtol, tol)
    atol = something(atol, tol)
    M = m_max + 2  # maximum m index needed
    # QuadGK requires an Array result container; Bumper output is exposed as a non-owning alias.
    result = Γₙ isa Array ? Γₙ : unsafe_wrap(Array, pointer(Γₙ), size(Γₙ))
    # Moment (y-d)^M peaks at |y-d| = √(M/2) and decays past it at least as fast as the Gaussian.
    ym = sqrt(T(M) / 2) + sqrt(-log(rtol))
    ymin = max(zero(T), d - ym)
    ymax = ym + d
    quadgk!(result, ymin, ymax; rtol, atol) do r, y
        ay = a * y
        jn = besselj(n, ay)
        g = iszero(ay) ? T(n == 1) / 2 : n * jn / ay
        V = (y * g, y * (besselj(n - 1, ay) - g), jn)
        prods = (V[1] * V[1], V[1] * V[2], V[1] * V[3], V[2] * V[2], V[2] * V[3], V[3] * V[3])
        w = exp(-(y - d)^2)
        yd = y - d
        @inbounds for m in 1:(M + 1)
            for q in 1:6
                r[2q - 1, m] = prods[q] * w
                r[2q, m] = prods[q] * w * y
            end
            w *= yd
        end
    end
    return Γₙ
end
