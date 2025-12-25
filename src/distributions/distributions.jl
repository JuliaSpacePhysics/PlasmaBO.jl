abstract type AbstractDistribution end

const PARTICLE_TYPE = Union{Symbol}

function _charge_mass(p, Z, A, q, m)
    # Resolve charge and mass based on particle type, Z, and A
    pp = CP.particle(p)
    q_val = @something(q, isnothing(Z) ? _charge(pp) : Z * qe)
    m_val = @something(m, isnothing(A) ? _mass(pp) : A * mp)
    return q_val, m_val
end

_mass(p::ParticleLike) = mass(p) / u"kg"
_charge(p::ParticleLike) = charge(p) / u"C"
_mass(s::AbstractDistribution) = s.m
_charge(s::AbstractDistribution) = s.q
_T(T::Quantity) = T / u"eV"
_T(T) = T
# charge(s) = s.q
# mass(s) = s.m
temperature(T) = T * q / kb
Debye_length(s::AbstractDistribution) = Debye_length(s.n, s.Tz * q / kb)
gyrofrequency(s) = gyrofrequency(_charge(s), s.B, _mass(s))
gyrofrequency(B, s) = gyrofrequency(_charge(s), B, _mass(s))

include("Maxwellian.jl")
include("BiKappa.jl")

for f in (:Maxwellian, :BiKappa, :BiKappa2)
    @eval $f(p::ParticleLike, args...; kw...) = $f(args...; particle = p, kw...)
end
