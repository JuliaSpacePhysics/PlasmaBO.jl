struct Maxwellian{T} <: AbstractDistribution
    q::T           # Charge (units of e)
    m::T           # Mass (units of mp)
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    vdr::T         # Perpendicular drift velocity (c)
end


Maxwellian(q, m, n, Tz, Tp = Tz; vdz = 0.0, vdr = 0.0) = Maxwellian(promote(q, m, n, Tz, Tp, vdz, vdr)...)
