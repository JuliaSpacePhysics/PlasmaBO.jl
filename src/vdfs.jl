"""
    HermiteVDF(; w_perp, w_para, alm, vd_para = 0, vr_perp = 0)

A distribution given by its Hermite–Hermite expansion coefficients `alm` about the
shifted-Gaussian reference with 1/e widths `w_perp, w_para` centered on `vd_para, vr_perp` [m/s],
see [`hermite_expansion`](@ref). [`BOHH`](@ref) uses it as is; every other model is expanded into one.
"""
struct HermiteVDF{T, A <: AbstractMatrix}
    w_perp::T
    w_para::T
    vd_para::T
    vr_perp::T
    alm::A
end

HermiteVDF(; w_perp, w_para, alm, vd_para = zero(w_para), vr_perp = zero(w_perp)) =
    HermiteVDF(promote(map(si_velocity, (w_perp, w_para, vd_para, vr_perp))...)..., alm)

# Tabulates `pdf(vdf, v⊥, v∥)` as `fv[i_para, i_perp]` on the expansion grid of `reference`.
function gen_fv2d(vdf, reference = hermite_reference(vdf))
    (; w_perp, w_para, vd_para, vr_perp) = reference
    v_para = range(vd_para - 10 * w_para, vd_para + 10 * w_para, step = 0.05 * w_para)
    v_perp = range(0, vr_perp + 10 * w_perp, step = 0.05 * w_perp)
    fv = [pdf(vdf, x, z) for z in v_para, x in v_perp]
    return (; v_para, v_perp, fv, w_perp, w_para, vd_para, vr_perp)
end

"""
    hermite_reference(vdf) -> (; w_perp, w_para, vd_para, vr_perp)

Shifted-Gaussian reference [m/s] that `vdf` is expanded about: its thermal speeds and mean drift,
so a model and its tabulation expand alike.
"""
hermite_reference(d::GaussianRing) = (; d.w_perp, w_para = d.vth_para, d.vd_para, d.vr_perp)
function hermite_reference(d::AbstractVDF)
    w_perp, w_para = vth_perp_para(d)
    vd_para = mean(d)[3]
    return (; w_perp, w_para, vd_para, vr_perp = zero(vd_para))
end

"""
    HermiteVDF(vdf::AbstractVDF; w_perp, w_para, vd_para, vr_perp, N_para = 16, N_perp = 16, atol = 0)

Expand a PlasmaPopulations model, about [`hermite_reference`](@ref) unless overridden.
"""
function HermiteVDF(
        vdf::AbstractVDF; w_perp = nothing, w_para = nothing, vd_para = nothing, vr_perp = nothing,
        N_para = 16, N_perp = 16, atol = 0.0,
    )
    isgyrotropic(vdf) || throw(ArgumentError("HermiteVDF expands gyrotropic models, not $(nameof(typeof(vdf)))"))
    ref = hermite_reference(vdf)
    g = gen_fv2d(
        vdf, (;
            w_perp = @something(w_perp, ref.w_perp), w_para = @something(w_para, ref.w_para),
            vd_para = @something(vd_para, ref.vd_para), vr_perp = @something(vr_perp, ref.vr_perp),
        )
    )
    alm = hermite_expansion(g.fv, g.v_para, g.v_perp, g.w_para, g.w_perp; N_para, N_perp, g.vd_para, g.vr_perp, atol).alm
    return HermiteVDF(g.w_perp, g.w_para, g.vd_para, g.vr_perp, alm)
end
# Exactly their own reference, the single term a₀₀ = 1. Keywords reach the numerical method above.
function HermiteVDF(d::Union{BiMaxwellian, GaussianRing})
    (; w_perp, w_para, vd_para, vr_perp) = hermite_reference(d)
    return HermiteVDF(w_perp, w_para, vd_para, vr_perp, ones(typeof(w_perp), 1, 1))
end
HermiteVDF(d::HermiteVDF) = d
HermiteVDF(d::Cold) = throw(ArgumentError("BOHH needs a thermal spread; use BOFluid for a cold population"))
