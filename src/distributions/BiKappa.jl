struct BiKappa{T <: Number} <: AbstractDistribution
    n::T
    κ::T
    Tz::T
    Tp::T
    vdz::T
    vdr::T
    q::T
    m::T
end

function BiKappa(n, κ, Tz, Tp = Tz; vdz = 0.0, vdr = 0.0, Z = nothing, A = nothing, q = nothing, m = nothing, particle = :p)
    q, m = _charge_mass(particle, Z, A, q, m)
    return BiKappa(promote(n, κ, Tz, Tp, vdz, vdr, q, m)...)
end

# BiKappa with different κ for z and x
struct BiKappa2{T <: Number} <: AbstractDistribution
    n::T
    κz::T
    κx::T
    Tz::T
    Tp::T
    vdz::T
    vdr::T
    q::T
    m::T
end

function BiKappa2(n, κz, κx, Tz, Tp = Tz; vdz = 0.0, vdr = 0.0, Z = nothing, A = nothing, q = nothing, m = nothing, particle = :p)
    q, m = _charge_mass(particle, Z, A, q, m)
    return BiKappa2(promote(n, κz, κx, Tz, Tp, vdz, vdr, q, m)...)
end

function _velocity_grid(vtx, vtz, vdx, vdz; dvx = nothing, dvz = nothing)
    dvx = @something dvx 0.05 * vtx
    dvz = @something dvz 0.05 * vtz
    vz_range = range(-10 * vtz + vdz, 10 * vtz + vdz, step = dvz)
    vx_range = range(0, 10 * vtx, step = dvx)

    vz = [z for z in vz_range, _ in vx_range]
    vx = [x for _ in vz_range, x in vx_range]

    return (; vz, vx, dvz, dvx)
end

# Thermal velocities with kappa correction
_vtp(s::BiKappa) = sqrt(2 * (1 - 1.5 / s.κ) * qe * s.Tp / s.m)
_vtz(s::BiKappa) = sqrt(2 * (1 - 1.5 / s.κ) * qe * s.Tz / s.m)

function gen_fv2d(s::BiKappa; dvx = nothing, dvz = nothing)
    κ = s.κ
    vdz = s.vdz * c0  # c -> m/s
    vdx = s.vdr * c0
    vtz = _vtz(s)
    vtx = _vtp(s)

    (; vz, vx, dvz, dvx) = _velocity_grid(vtx, vtz, vdx, vdz; dvx, dvz)

    coef = 1 / (sqrt(π^3 * κ^3) * vtz * vtx^2) * gamma(κ + 1) / gamma(κ - 0.5)
    fv = @. coef * (1 + (vz - vdz)^2 / (κ * vtz^2) + vx^2 / (κ * vtx^2))^(-κ - 1)

    return (; vz, vx, fv, dvz, dvx, vtz, vtx, vdz, vdx)
end

_vtz(s::BiKappa2) = sqrt(2 * (1 - 0.5 / s.κz) * q * s.Tz / _mass(s))
_vtp(s::BiKappa2) = sqrt(2 * (1 - 1 / s.κx) * q * s.Tp / _mass(s))

function gen_fv2d(s::BiKappa2; dvx = nothing, dvz = nothing)
    κz = s.κz
    κx = s.κx
    vdz = s.vdz * c0  # c -> m/s
    vdx = s.vdr * c0
    vtz = _vtz(s)
    vtx = _vtp(s)

    (; vz, vx, dvz, dvx) = _velocity_grid(vtx, vtz, vdx, vdz; dvx, dvz)

    coef = 1 / (sqrt(π^3 * κz) * vtz * vtx^2) * exp(loggamma(κz + 1) - loggamma(κz + 0.5))
    fv = @. coef * (1 + (vz - vdz)^2 / (κz * vtz^2))^(-κz - 1) * (1 + vx^2 / (κx * vtx^2))^(-κx - 1)

    return (; vz, vx, fv, dvz, dvx, vtz, vtx, vdz, vdx)
end
