"""
    BOHH(; N = 2, J = 8)

Dispersion solver using the Hermite-Hermite (HH) basis expansion.

`N` controls the truncation order of the cyclotron harmonic index.
`J` controls the truncation order of the number of poles for Z-function approximation

This method transforms the dispersion relation into a matrix eigenvalue problem
using J-pole approximation for the plasma dispersion function, allowing
simultaneous computation of all wave modes.
"""
@kwdef struct BOHH <: AbstractDispersionAlgorithm
    N::Int = 2
    J::Int = 8
end

struct HHSolverParam{T}
    wc::T               # signed cyclotron frequency
    wp::T
    vtz::T              # parallel, perpendicular reference widths
    vtp::T
    vdz::T              # parallel drift
    vdr::T              # perpendicular ring speed
    alm::Matrix{T}      # power-basis Hermite coefficients
end

# Susceptibility is linear in alm: a negative density is |n| with negated coefficients.
function HHSolverParam(s, ::Type{T}) where {T}
    iszero(s.Omega) && throw(ArgumentError("BOHH expands about the gyration: needs B0 ≠ 0"))
    h = HermiteVDF(s.vdf)
    c = T(s.c)
    return HHSolverParam{T}(
        T(s.Omega), sqrt(abs(T(s.Pi2))), T(h.w_para) / c, T(h.w_perp) / c, T(h.vd_para) / c, T(h.vr_perp) / c,
        sign(T(s.Pi2)) * Matrix{T}(h.alm)
    )
end

# alm-weighted m-sums of the `_compute_integral_matrices!` moments, for each l and product q:
#   Tₗ[q, l+1] = Σₘ alm[l+1,m+1] · (2·Kq(m+1) - m·Kq(m-1))   (Kq = Γ[2q-1, :], the lower power of y)
#   Uₗ[q, l+1] = Σₘ alm[l+1,m+1] · Γ[2q, m+1]
function _precompute_lm_sums!(Tₗ, Uₗ, alm, Γ)
    l_max, m_max = size(alm) .- 1
    @inbounds for l in 0:l_max, q in 1:6
        a = alm[l + 1, 1]
        t = a * 2 * Γ[2q - 1, 2]
        u = a * Γ[2q, 1]
        for m in 1:m_max
            a = alm[l + 1, m + 1]
            t += a * (2 * Γ[2q - 1, m + 2] - m * Γ[2q - 1, m])
            u += a * Γ[2q, m + 1]
        end
        Tₗ[q, l + 1] = t
        Uₗ[q, l + 1] = u
    end
    return nothing
end

# Per product q: P = Σₗ cˡ Tₗ[q, l+1] and Q = Σₗ (2cˡ⁺¹ - l·cˡ⁻¹) Uₗ[q, l+1] = 2c·R(c) - R′(c),
# with R(c) = Σₗ cˡ Uₗ[q, l+1], by Horner.
@inline function _pole_sums(Tₗ, Uₗ, c)
    z = zero(c)
    P = R = dR = ntuple(_ -> z, 6)
    @inbounds for l in size(Tₗ, 2):-1:1
        dR = dR .* c .+ R
        R = R .* c .+ ntuple(q -> Uₗ[q, l], 6)
        P = P .* c .+ ntuple(q -> Tₗ[q, l], 6)
    end
    return P, 2c .* R .- dR
end

# Fills poles ξ and E-couplings B at r0 + (n + N)·J + j for one species; returns its pole-free coupling.
function _assemble_species!(ξ, B, p::HHSolverParam{T}, r0, kx, kz, (; bzj, czj), N) where {T}
    (; wc, vtz, vtp, vdz, alm) = p
    J = length(czj)
    wp2 = p.wp^2
    as = kx * vtp / wc
    d = p.vdr / vtp
    R = exp(-d^2) + sqrt(T(π)) * d * erfc(-d)
    dr = vdz / vtz
    vq = vtp^2 / vtz
    l_max, m_max = size(alm) .- 1
    b = ntuple(_ -> zero(eltype(ξ)), 6)

    @no_escape begin
        Γ = @temp_array(T, 12, m_max + 3)
        Ils = @temp_array(T, l_max + 3)
        Tₗ = @temp_array(T, 6, l_max + 1)
        Uₗ = @temp_array(T, 6, l_max + 1)
        Ils .= funIn.(T, 0:(l_max + 2))

        # Γ serves both ±a: n enters otherwise only through n·wc and the sign of the V₁ products.
        for a in 0:N
            _compute_integral_matrices!(Γ, a, as, d, m_max)
            Γ .*= 2 / R
            _precompute_lm_sums!(Tₗ, Uₗ, alm, Γ)

            s11 = s22 = s33 = zero(T)
            @inbounds for l in 0:l_max
                Iₗ₋₁ = l >= 1 ? Ils[l] : zero(T)
                s11 += Ils[l + 1] * Tₗ[1, l + 1]
                s22 += Ils[l + 1] * Tₗ[4, l + 1]
                s33 += (dr * (2 * Ils[l + 2] - l * Iₗ₋₁) + (2 * Ils[l + 3] - l * Ils[l + 1])) * Uₗ[6, l + 1]
            end
            mult = a == 0 ? 1 : 2
            b = b .- mult * wp2 .* (s11, 0, 0, s22, 0, s33)

            for j in 1:J
                cⱼ = kz < 0 ? -czj[j] : czj[j]   # the poles of Z(ζ) mirror for k∥ < 0
                P, Q = _pole_sums(Tₗ, Uₗ, cⱼ)
                u = cⱼ * vtz + vdz
                w = u / vtp
                for n in (a == 0 ? (0,) : (-a, a))
                    cnj = n * wc + kz * u
                    # At n = 0 the pole's k∥ cancels exactly, keeping k∥ = 0 regular.
                    X = (wp2 * bzj[j]) .* (iszero(n) ? (vq / u) .* Q : (n * wc .* P .+ (kz * vq) .* Q) ./ cnj)
                    sn = sign(n)
                    pj = (X[1], im * sn * X[2], sn * w * X[3], X[4], -im * w * X[5], w^2 * X[6])
                    r = r0 + (n + N) * J + j
                    ξ[r] = cnj
                    _coupling!(B, r, pj)
                    b = b .- pj
                end
            end
        end
    end
    return b
end

# Velocity triples (x, y, z) at pole ξ couple to E through
# [p11 p12 p13; -p12 p22 p23; p13 -p23 p33], and E to them through -I (the current).
function _coupling!(B, r, (p11, p12, p13, p22, p23, p33))
    B[1, 1, r], B[1, 2, r], B[1, 3, r] = p11, p12, p13
    B[2, 1, r], B[2, 2, r], B[2, 3, r] = -p12, p22, p23
    B[3, 1, r], B[3, 2, r], B[3, 3, r] = p13, -p23, p33
    return B
end

prepare(alg::BOHH, populations, ::Type{T}) where {T} =
    (; species = map(sp -> HHSolverParam(sp, T), populations), jpole = get_jpole_coefficients(alg.J, T))

# One triple per (species, n, j) pole, then one pole-free triple (ξ = 0) per species.
function _operator(alloc, params, k::Wavevector, alg::BOHH, ::Type{CT}) where {CT}
    (; species, jpole) = params
    S, J, N = length(species), length(jpole.czj), alg.N
    SNJ = S * (2N + 1) * J
    m = SNJ + S
    ξ = alloc(CT, m)
    B = alloc(CT, 3, 3, m)
    tol = alloc(real(CT), m)
    kx, kz = real(CT)(k.k_perp), real(CT)(k.k_para)
    for (s, param) in enumerate(species)
        r0 = (s - 1) * (2N + 1) * J
        b = _assemble_species!(ξ, B, param, r0, kx, kz, jpole, N)
        ξ[SNJ + s] = 0
        _coupling!(B, SNJ + s, b)
        # Relative to each species' largest residue: a global one would let electrons erase ion harmonics.
        rs = (r0 + 1):(r0 + (2N + 1) * J)
        t = 3eps(real(CT)) * max(norm(view(B, :, :, SNJ + s)), maximum(r -> norm(view(B, :, :, r)), rs))
        tol[rs] .= t
        tol[SNJ + s] = t
    end
    return _reduced_operator(alloc, ξ, B, tol, kx, kz)
end
