module PlasmaBOSpecialFunctionsExt

using PlasmaBO
using SpecialFunctions

PlasmaBO.besselj(n::Integer, x) = SpecialFunctions.besselj(n, x)
PlasmaBO._erfc(x) = SpecialFunctions.erfc(x)
PlasmaBO.loggamma(x::Number) = SpecialFunctions.loggamma(x)
PlasmaBO.gamma(x::Number) = SpecialFunctions.gamma(x)

end
