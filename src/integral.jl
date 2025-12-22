"""
    funAn(n, a, d, m, p)

Perpendicular integral for J_n²: ∫ J_n²(ay) exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funAn(n, a, d, m, p; rtol = 1.0e-10, atol = 1.0e-10)
    m < 0 && return 0.0

    ym = 10.0
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

    ym = 10.0
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

    ym = 10.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    function integrand(y)
        jn_deriv = 0.5 * (besselj(n - 1, a * y) - besselj(n + 1, a * y))
        return 0.25 * jn_deriv^2 * exp(-(y - d)^2) * (y - d)^m * y^p
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
