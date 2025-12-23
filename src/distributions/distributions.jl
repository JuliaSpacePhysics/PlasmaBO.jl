abstract type AbstractDistribution end

charge(s) = s.q * q
mass(s) = s.m * mp
temperature(T) = T * q / kb
Debye_length(s::AbstractDistribution) = Debye_length(s.n, s.Tz * q / kb)

include("Maxwellian.jl")
