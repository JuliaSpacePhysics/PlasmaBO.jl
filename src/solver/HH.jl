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

# alm-weighted m-sums of the `_compute_integral_matrices!` moments, for each l and k ∈ (A, B, C):
#   Tₗ[k, l+1] = Σₘ alm[l+1,m+1] · (2·Kₖ(m+1) - m·Kₖ(m-1))   (Kₖ = Γ[2k-1, :], the lower power of y)
#   Uₗ[k, l+1] = Σₘ alm[l+1,m+1] · Γ[2k, m+1]
function _precompute_lm_sums!(Tₗ, Uₗ, alm, Γ)
    l_max, m_max = size(alm) .- 1
    @inbounds for l in 0:l_max
        a = alm[l + 1, 1]
        ta = a * 2 * Γ[1, 2]
        ua = a * Γ[2, 1]
        tb = a * 2 * Γ[3, 2]
        ub = a * Γ[4, 1]
        tc = a * 2 * Γ[5, 2]
        uc = a * Γ[6, 1]
        for m in 1:m_max
            a = alm[l + 1, m + 1]
            ta += a * (2 * Γ[1, m + 2] - m * Γ[1, m])
            ua += a * Γ[2, m + 1]
            tb += a * (2 * Γ[3, m + 2] - m * Γ[3, m])
            ub += a * Γ[4, m + 1]
            tc += a * (2 * Γ[5, m + 2] - m * Γ[5, m])
            uc += a * Γ[6, m + 1]
        end
        Tₗ[1, l + 1] = ta; Tₗ[2, l + 1] = tb; Tₗ[3, l + 1] = tc
        Uₗ[1, l + 1] = ua; Uₗ[2, l + 1] = ub; Uₗ[3, l + 1] = uc
    end
    return nothing
end

# Per k: P = Σₗ cˡ Tₗ[k, l+1] and Q = Σₗ (2cˡ⁺¹ - l·cˡ⁻¹) Uₗ[k, l+1] = 2c·R(c) - R′(c),
# with R(c) = Σₗ cˡ Uₗ[k, l+1], by Horner.
@inline function _pole_sums(Tₗ, Uₗ, c)
    z = zero(c)
    P = R = dR = (z, z, z)
    @inbounds for l in size(Tₗ, 2):-1:1
        dR = dR .* c .+ R
        R = R .* c .+ (Uₗ[1, l], Uₗ[2, l], Uₗ[3, l])
        P = P .* c .+ (Tₗ[1, l], Tₗ[2, l], Tₗ[3, l])
    end
    return P, 2c .* R .- dR
end

# Rows r, r + SNJ1, r + 2SNJ1 couple to E through [p11 p12 p13; -p12 p22 p23; p13 -p23 p33].
function _couple_E!(M, r, SNJ1, SNJ3, (p11, p12, p13, p22, p23, p33))
    M[r, (SNJ3 + 1):(SNJ3 + 3)] .= (p11, p12, p13)
    M[r + SNJ1, (SNJ3 + 1):(SNJ3 + 3)] .= (-p12, p22, p23)
    M[r + 2SNJ1, (SNJ3 + 1):(SNJ3 + 3)] .= (p13, -p23, p33)
    return M
end

# Fills the pole rows r0 + (n + N)·J + j of one species; returns its E-coupling for the species rows.
function _assemble_species!(M, p::HHSolverParam{T}, r0, kx, kz, SNJ1, SNJ3, czj, bzj, N) where {T}
    (; wc, vtz, vtp, vdz, alm) = p
    J = length(czj)
    wp2 = p.wp^2
    as = kx * vtp / wc
    d = p.vdr / vtp
    R = exp(-d^2) + sqrt(T(π)) * d * erfc(-d)
    dr = vdz / vtz
    vq = vtp^2 / vtz
    l_max, m_max = size(alm) .- 1
    b = ntuple(_ -> zero(eltype(M)), 6)

    @no_escape begin
        Γ = @temp_array(T, 6, m_max + 3)
        Ils = @temp_array(T, l_max + 3)
        Tₗ = @temp_array(T, 3, l_max + 1)
        Uₗ = @temp_array(T, 3, l_max + 1)
        Ils .= funIn.(T, 0:(l_max + 2))

        # Everything per harmonic depends on n only through n·wc and Γ, which is even in n.
        for a in 0:N
            _compute_integral_matrices!(Γ, a, as, d, m_max)
            Γ .*= 2 / R
            _precompute_lm_sums!(Tₗ, Uₗ, alm, Γ)

            s11 = s22 = s33 = zero(T)
            @inbounds for l in 0:l_max
                Iₗ₋₁ = l >= 1 ? Ils[l] : zero(T)
                s11 += Ils[l + 1] * Tₗ[1, l + 1]
                s22 += Ils[l + 1] * Tₗ[3, l + 1]
                s33 += (dr * (2 * Ils[l + 2] - l * Iₗ₋₁) + (2 * Ils[l + 3] - l * Ils[l + 1])) * Uₗ[1, l + 1]
            end
            nwkp = a * wc / (kx * vtp)
            mult = a == 0 ? 1 : 2
            b = b .- mult * wp2 .* (nwkp^2 * s11, 0, 0, s22, 0, s33)

            for j in 1:J
                cⱼ = czj[j]
                P, Q = _pole_sums(Tₗ, Uₗ, cⱼ)
                u = cⱼ * vtz + vdz
                w = u / vtp
                for n in (a == 0 ? (0,) : (-a, a))
                    cnj = n * wc + kz * u
                    # At n = 0 the pole's k∥ cancels exactly, keeping k∥ = 0 regular.
                    X = (wp2 * bzj[j]) .* (iszero(n) ? (vq / u) .* Q : (n * wc .* P .+ (kz * vq) .* Q) ./ cnj)
                    nwkp = n * wc / (kx * vtp)
                    pj = (nwkp^2 * X[1], im * nwkp * X[2], nwkp * w * X[1], X[3], -im * w * X[2], w^2 * X[1])
                    r = r0 + (n + N) * J + j
                    M[r, r] = M[r + SNJ1, r + SNJ1] = M[r + 2SNJ1, r + 2SNJ1] = cnj
                    _couple_E!(M, r, SNJ1, SNJ3, pj)
                    b = b .- pj
                end
            end
        end
    end
    return b
end

# State Vector Organization
# Indices 1 to SNJ1:           v_snj_x components + species j_x
# Indices SNJ1+1 to 2*SNJ1:    v_snj_y components + species j_y
# Indices 2*SNJ1+1 to 3*SNJ1:  v_snj_z components + species j_z
# Indices SNJ3+1 to SNJ3+6:    E_x, E_y, E_z, B_x, B_y, B_z
# where SNJ1 = SNJ + S (SNJ pole velocities + S species auxiliary j's)

matrix_size(alg::BOHH, species) = let S = length(species)
    3 * (S * (2 * alg.N + 1) * alg.J + S) + 6
end

prepare(::BOHH, populations, ::Type{T}) where {T} = map(sp -> HHSolverParam(sp, T), populations)

function dispersion_matrix!(M, params, k::Wavevector, alg::BOHH)
    T = real(eltype(M))
    N = alg.N
    S = length(params)
    (; J, bzj, czj) = get_jpole_coefficients(alg.J, T)

    # The k⊥ → 0 limit is finite but each factor (nwkp, Bessel moments) is singular there.
    kx = iszero(k.k_perp) ? T(1.0e-30) : T(k.k_perp)
    kz = T(k.k_para)
    kz < 0 && (czj = -czj)
    SNJ = S * (2N + 1) * J
    SNJ1 = SNJ + S
    SNJ3 = 3SNJ1

    for (s, param) in enumerate(params)
        b = _assemble_species!(M, param, (s - 1) * (2N + 1) * J, kx, kz, SNJ1, SNJ3, czj, bzj, N)
        _couple_E!(M, SNJ + s, SNJ1, SNJ3, b)
    end

    # E(J) coupling: J_xyz = j_xyz + sum(v_snj_xyz)
    M[SNJ3 + 1, 1:SNJ1] .= -1
    M[SNJ3 + 2, (SNJ1 + 1):2SNJ1] .= -1
    M[SNJ3 + 3, (2SNJ1 + 1):3SNJ1] .= -1
    _B_E_part!(M, SNJ3, kx, kz)
    return M
end

# Maxwell's curl equations in units with c = 1
function _B_E_part!(M, idx, kx, kz)
    # E(B) coupling: Ampère's law
    M[idx + 1, idx + 5] = kz
    M[idx + 2, idx + 4] = -kz
    M[idx + 2, idx + 6] = kx
    M[idx + 3, idx + 5] = -kx

    # B(E) coupling: Faraday's law
    M[idx + 4, idx + 2] = -kz
    M[idx + 5, idx + 1] = kz
    M[idx + 5, idx + 3] = -kx
    M[idx + 6, idx + 2] = kx
    return M
end
