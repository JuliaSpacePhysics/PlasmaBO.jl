using Bessels: besselj

"""
    funIn(n)

Normalization integral I_n = ∫ v^n exp(-v²) dv / √π.

Returns Γ((n+1)/2)/√π for even n, 0 for odd n.
"""
function funIn(n)
    return isodd(n) ? 0.0 : gamma(0.5 * (n + 1)) / sqrt(π)
end

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
# Layout — `Γₙ` is `(6, m_max+4)`; column 1 is unused padding (so callers can read
# `Γₙ[k, m+1]` for m=0 without bounds gymnastics) and `Γₙ[:, m+2]` holds the six
# integrals at moment `m` (m = 0:m_max+2):
#   Γₙ[1, m+2] = Aₙ(m, 0)   Γₙ[2, m+2] = Aₙ(m, 1)
#   Γₙ[3, m+2] = Bₙ(m, 1)   Γₙ[4, m+2] = Bₙ(m, 2)
#   Γₙ[5, m+2] = Cₙ(m, 2)   Γₙ[6, m+2] = Cₙ(m, 3)
function _compute_integral_matrices!(Γₙ, n, a, d, m_max; rtol = 1.0e-10, atol = 1.0e-10)
    M = m_max + 2  # maximum m index needed
    # View columns 2:end as a contiguous Array for quadgk! to fill directly,
    # skipping the padding column 1.
    vals = unsafe_wrap(Array, pointer(Γₙ) + sizeof(eltype(Γₙ)) * 6, (6, M + 1))
    ym = 5.0
    ymin = max(0.0, d - ym)
    ymax = ym + d
    quadgk!(vals, ymin, ymax; rtol, atol) do r, y
        ay = a * y
        jn = besselj(n, ay)
        jnm = besselj(n - 1, ay)
        # J_{n+1} from forward recurrence; integrand contribution at ay≈0 is
        # negligible so the 1-bit cancellation there does not affect the integral.
        jnp = iszero(ay) ? besselj(n + 1, ay) : (2 * n / ay) * jn - jnm
        jnd = 0.5 * (jnm - jnp)

        jn2 = jn * jn
        jnjd = jn * jnd
        jnd2 = jnd * jnd

        e = exp(-(y - d)^2)
        yd = y - d
        y2 = y * y
        y3 = y2 * y

        ydm = 1.0
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
