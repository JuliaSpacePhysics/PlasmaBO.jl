using Bessels: Bessels
using Gamma: Gamma

# Normalization integral Iₙ = ∫ vⁿ exp(-v²) dv / √π: Γ((n+1)/2)/√π for even n, 0 for odd n.
function funIn(::Type{T}, n) where {T <: Real}
    isodd(n) && return zero(T)
    I = one(T)
    for k in 2:2:n
        I *= T(k - 1) / T(2)
    end
    return I
end

# `SpecialFunctions` extensions override for `BigFloat` and `Double64`
besselj(n::Integer, x::Union{Float32, Float64}) = Bessels.besselj(n, x)
# Complementary error function via libopenlibm
# Same path as SpecialFunctions for Float64
erfc(x::Number) = _erfc(float(x))

_erfc(x::Float64) = ccall((:erfc, Base.Math.libm), Float64, (Float64,), x)
_erfc(x::Float32) = ccall((:erfcf, Base.Math.libm), Float32, (Float32,), x)

# Gamma.jl covers Float16/32/64, BigFloat, Integer, Rational.
# SpecialFunctions extends coverage to anything else (e.g. Double64)
const _GammaSupported = Union{Float16, Float32, Float64, BigFloat, Integer, Rational}
loggamma(x::_GammaSupported) = Gamma.loggamma(x)
gamma(x::_GammaSupported) = Gamma.gamma(x)
