using Bessels: besselj

"""
    funAn(n, a, d, m, p)

Perpendicular integral for J_n²: ∫ J_n²(ay) exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funAn(n, a, d, m, p; rtol = 1.0e-10, atol = 1.0e-10)
    m < 0 && return 0.0

    ym = 5.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    integrand(y) = besselj(n, a * y)^2 * exp(-(y - d)^2) * (y - d)^m * y^p
    return quadgk(integrand, ymin, ymax; rtol, atol)[1]
end


"""
    funBn(n, a, d, m, p)

Perpendicular integral for J_n·J_n': ∫ J_n(ay) J_n'(ay) exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funBn(n, a, d, m, p; rtol = 1.0e-10, atol = 1.0e-10)
    m < 0 && return 0.0

    ym = 5.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    function integrand(y)
        jn = besselj(n, a * y)
        jn_deriv = 0.5 * (besselj(n - 1, a * y) - besselj(n + 1, a * y))
        return jn * jn_deriv * exp(-(y - d)^2) * (y - d)^m * y^p
    end

    return quadgk(integrand, ymin, ymax; rtol, atol)[1]
end

"""
    funCn(n, a, d, m, p)

Perpendicular integral for (J_n')²: ∫ [J_n'(ay)]² exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funCn(n, a, d, m, p; rtol = 1.0e-10, atol = 1.0e-10)
    m < 0 && return 0.0

    ym = 5.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    function integrand(y)
        jn_deriv = 0.5 * (besselj(n - 1, a * y) - besselj(n + 1, a * y))
        return jn_deriv^2 * exp(-(y - d)^2) * (y - d)^m * y^p
    end

    return quadgk(integrand, ymin, ymax; rtol, atol)[1]
end


"""
    funIn(n)

Normalization integral I_n = ∫ v^n exp(-v²) dv / √π.

Returns Γ((n+1)/2)/√π for even n, 0 for odd n.
"""
function funIn(n)
    return isodd(n) ? 0.0 : gamma(0.5 * (n + 1)) / sqrt(π)
end

# Fill pre-allocated matrices `Aₙ`, `Bₙ`, `Cₙ` for harmonic index `n` and all m in
# 0:(m_max+2) using a **single** adaptive quadrature call.
# This is equivalent to calling `funAn`, `funBn`, `funCn` separately for every m but
# evaluates the three Bessel functions J_{n-1}, J_n, J_{n+1} only once per quadrature
# node, reducing the number of `quadgk` calls from 6*(m_max+3) to 1.

# Fills:
#   Aₙ[m+2, 1] = funAn(n, a, d, m, 0),  Aₙ[m+2, 2] = funAn(n, a, d, m, 1)
#   Bₙ[m+2, 1] = funBn(n, a, d, m, 1),  Bₙ[m+2, 2] = funBn(n, a, d, m, 2)
#   Cₙ[m+2, 1] = funCn(n, a, d, m, 2),  Cₙ[m+2, 2] = funCn(n, a, d, m, 3)
function _compute_integral_matrices!(Aₙ, Bₙ, Cₙ, n, a, d, m_max; kw...)
    M = m_max + 2  # maximum m index needed
    return @no_escape begin
        _vals_raw = @alloc(Float64, 6, M + 1)
        vals = unsafe_wrap(Array, pointer(_vals_raw), (6, M + 1))
        _compute_integral_matrices!(vals, n, a, d, m_max; kw...)
        @inbounds for m in 0:M
            base = 6 * m
            idm = m + 2
            Aₙ[idm, 1] = vals[base + 1]
            Aₙ[idm, 2] = vals[base + 2]
            Bₙ[idm, 1] = vals[base + 3]
            Bₙ[idm, 2] = vals[base + 4]
            Cₙ[idm, 1] = vals[base + 5]
            Cₙ[idm, 2] = vals[base + 6]
        end
    end
end

function _compute_integral_matrices!(vals, n, a, d, m_max; rtol = 1.0e-10, atol = 1.0e-10)
    ym = 5.0
    ymin = max(0.0, d - ym)
    ymax = ym + d
    M = m_max + 2  # maximum m index needed
    return quadgk!(vals, ymin, ymax; rtol, atol) do r, y
        ay = a * y
        jn = besselj(n, ay)
        jnm = besselj(n - 1, ay)
        jnp = besselj(n + 1, ay)
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
            r[1, m] = jn2 * e * ydm        # funAn(m, 0)
            r[2, m] = jn2 * e * ydm * y    # funAn(m, 1)
            r[3, m] = jnjd * e * ydm * y    # funBn(m, 1)
            r[4, m] = jnjd * e * ydm * y2   # funBn(m, 2)
            r[5, m] = jnd2 * e * ydm * y2   # funCn(m, 2)
            r[6, m] = jnd2 * e * ydm * y3   # funCn(m, 3)
            ydm *= yd
        end
    end

end
