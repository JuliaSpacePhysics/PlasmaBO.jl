# Fill `Γₙ` for harmonic index `n` and all m in 0:(m_max+2) with a single adaptive
# quadrature call. Evaluates J_{n-1}, J_n only once per quadrature node and derives
# J_{n+1} by recurrence, replacing 6·(m_max+3) separate `quadgk` calls with 1.
#
# The six integrals (all over y ∈ [max(0, d-5), d+5] with weight exp(-(y-d)²)):
#   Aₙ(m, p) = ∫ J_n(ay)²              · exp(-(y-d)²) · (y-d)^m · y^p dy
#   Bₙ(m, p) = ∫ J_n(ay) · J_n'(ay)    · exp(-(y-d)²) · (y-d)^m · y^p dy
#   Cₙ(m, p) = ∫ J_n'(ay)²             · exp(-(y-d)²) · (y-d)^m · y^p dy
# with J_n'(x) = (J_{n-1}(x) - J_{n+1}(x)) / 2.
#
# Layout — `Γₙ` is `(6, m_max+3)` and `Γₙ[:, m+1]` holds the six
# integrals at moment `m` (m = 0:m_max+2):
#   Γₙ[1, m+1] = Aₙ(m, 0)   Γₙ[2, m+1] = Aₙ(m, 1)
#   Γₙ[3, m+1] = Bₙ(m, 1)   Γₙ[4, m+1] = Bₙ(m, 2)
#   Γₙ[5, m+1] = Cₙ(m, 2)   Γₙ[6, m+1] = Cₙ(m, 3)
function _compute_integral_matrices!(Γₙ, n, a, d, m_max; rtol = nothing, atol = nothing)
    T = eltype(Γₙ)
    tol = T === Float64 ? T(1.0e-10) : sqrt(eps(one(T)))
    rtol = something(rtol, tol)
    atol = something(atol, tol)
    M = m_max + 2  # maximum m index needed
    # QuadGK requires an Array result container. Managed `@temp_array` output
    # is already an Array; Bumper output is exposed as a non-owning Array alias.
    result = Γₙ isa Array ? Γₙ : unsafe_wrap(Array, pointer(Γₙ), size(Γₙ))
    ym = T(5)
    ymin = max(zero(T), d - ym)
    ymax = ym + d
    quadgk!(result, ymin, ymax; rtol, atol) do r, y
        ay = a * y
        jn = besselj(n, ay)
        jnm = besselj(n - 1, ay)
        # J_{n+1} from forward recurrence; integrand contribution at ay≈0 is
        # negligible so the 1-bit cancellation there does not affect the integral.
        jnp = iszero(ay) ? besselj(n + 1, ay) : (T(2) * n / ay) * jn - jnm
        jnd = (jnm - jnp) / T(2)

        jn2 = jn * jn
        jnjd = jn * jnd
        jnd2 = jnd * jnd

        e = exp(-(y - d)^2)
        yd = y - d
        y2 = y * y
        y3 = y2 * y

        ydm = one(T)
        @inbounds for m in 1:(M + 1)
            r[1, m] = jn2 * e * ydm        # Aₙ(m-1, 0)
            r[2, m] = jn2 * e * ydm * y    # Aₙ(m-1, 1)
            r[3, m] = jnjd * e * ydm * y   # Bₙ(m-1, 1)
            r[4, m] = jnjd * e * ydm * y2  # Bₙ(m-1, 2)
            r[5, m] = jnd2 * e * ydm * y2  # Cₙ(m-1, 2)
            r[6, m] = jnd2 * e * ydm * y3  # Cₙ(m-1, 3)
            ydm *= yd
        end
    end
    return Γₙ
end
