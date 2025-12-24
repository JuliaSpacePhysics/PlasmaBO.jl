struct Maxwellian{T <: Number} <: AbstractDistribution
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    vdr::T         # Perpendicular drift velocity (c)
    q::T           # Charge
    m::T           # Mass (units of mp)
end

Base.eltype(::Maxwellian{T}) where {T} = T

function Maxwellian(n, Tz, Tp = Tz; vdz = 0.0, vdr = 0.0, Z = nothing, A = nothing, q = nothing, m = nothing, particle = :p)
    q, m = _charge_mass(particle, Z, A, q, m)
    return Maxwellian(promote(n, _T(Tz), _T(Tp), vdz, vdr, q, m)...)
end


Maxwellian(p::ParticleLike, args...; kw...) = Maxwellian(args...; particle = p, kw...)
