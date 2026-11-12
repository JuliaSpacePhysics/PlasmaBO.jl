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
        Γ = _temp_array(T, 12, m_max + 3)
        Ils = _temp_array(T, l_max + 3)
        Tₗ = _temp_array(T, 6, l_max + 1)
        Uₗ = _temp_array(T, 6, l_max + 1)
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
                    _coupling!(view(B, :, :, r), pj)
                    b = b .- pj
                end
            end
        end
    end
    return b
end

prepare(alg::BOHH, populations, ::Type{T}) where {T} =
    (; species = map(sp -> HHSolverParam(sp, T), populations), jpole = get_jpole_coefficients(alg.J, T))

# One triple per (species, n, j) pole, then one pole-free triple (ξ = 0) per species; E couples
# to each through -I (the current).
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
        _coupling!(view(B, :, :, SNJ + s), b)
        # Relative to each species' largest residue: a global one would let electrons erase ion harmonics.
        rs = (r0 + 1):(r0 + (2N + 1) * J)
        t = 3eps(real(CT)) * max(norm(view(B, :, :, SNJ + s)), maximum(r -> norm(view(B, :, :, r)), rs))
        tol[rs] .= t
        tol[SNJ + s] = t
    end
    return _reduced_operator(alloc, ξ, B, tol, kx, kz)
end

# The operator with auxiliary triples r (ẋᵣ = ξᵣ xᵣ + Bᵣ E, Ė = -Σᵣ xᵣ + k × B) reduced to the states
# the fields couple to: at each distinct pole, as many as the rank of the summed residue, with pivots
# below the poles' `tol` cut (Tian, Xie & Wei 2026, arXiv:2610.41450, §3.1). Roots are unchanged
# (their eq. 15); the dropped eigenvalues sit exactly on the poles.
function _reduced_operator(alloc, ξ, B, tol, kx, kz)
    CT = eltype(ξ)
    m = length(ξ)
    # + 0 maps -0.0 to 0.0, so equal poles sort into contiguous runs.
    perm = sortperm(ξ; by = z -> (real(z) + 0, imag(z) + 0))
    ξr = alloc(CT, 3m)
    L = alloc(CT, 3, 3m)    # E rows ← reduced states
    Rf = alloc(CT, 3m, 3)   # reduced states ← E
    Q = alloc(CT, 3, 3)
    r = 0
    i = 1
    while i <= m
        pole = ξ[perm[i]]
        fill!(Q, zero(CT))
        cut = zero(real(CT))
        while i <= m && ξ[perm[i]] == pole
            Q .-= view(B, :, :, perm[i])
            cut = max(cut, tol[perm[i]])
            i += 1
        end
        # Complete-pivoting elimination of Q: each pivot Q[a, b] peels the rank-1 term Q[:, b] Q[a, :] / Q[a, b].
        for _ in 1:3
            p, idx = findmax(abs, Q)
            p > cut || break
            a, b = Tuple(idx)
            s = sqrt(Q[a, b])   # split the pivot's size evenly between the two couplings
            r += 1
            ξr[r] = pole
            for c in 1:3
                L[c, r] = Q[c, b] / s
                Rf[r, c] = Q[a, c] / s
            end
            for c in 1:3, d in 1:3
                Q[d, c] -= L[d, r] * Rf[r, c]
            end
        end
    end
    M = zero!(alloc(CT, r + 6, r + 6))
    for t in 1:r
        M[t, t] = ξr[t]
    end
    M[1:r, (r + 1):(r + 3)] .= view(Rf, 1:r, :)
    M[(r + 1):(r + 3), 1:r] .= view(L, :, 1:r)
    return _B_E_part!(M, r, kx, kz)
end

# D(ω) and D′(ω), D the E rows of the operator with its states and B eliminated:
# D(ω) = ωI - Σₜ Lₜ Rₜ/(ω - ξₜ) - M_EB M_BE/ω.
function _D_dD(M, r, K, ω)
    D = ω * I - K / ω
    dD = I + K / ω^2
    for t in 1:r
        d = ω - M[t, t]
        for b in 1:3, a in 1:3
            q = M[r + a, t] * M[t, r + b] / d
            D[a, b] -= q
            dD[a, b] += q / d
        end
    end
    return D, dD
end

# Newton on det D(ω) = 0, step δ = 1/tr(D⁻¹D′). D's terms are rounded one by one, so its roots keep the
# digits the dense eigensolver loses at large J. det D cancels, so convergence is judged on the step: below
# atol, or no longer halving below √atol. A full step from a poor seed can land in a neighbouring root's
# basin, so steps are damped by λ = min(1, 1/2h), h the last step's relative departure from its linear
# prediction (1 - λ)δ; λ = 1/2 at the seed, where h is unknown.
function _newton_det(M, r, K, ω, atol, maxiter)
    prev = oftype(abs(ω), Inf)
    λ = oftype(prev, 0.5)
    δ₀ = zero(ω)
    for it in 1:maxiter
        D, dD = _D_dD(M, r, K, ω)
        all(isfinite, D) && all(isfinite, dD) || return oftype(ω, NaN), it, ReturnCode.Failure
        F = lu!(D; check = false)
        issuccess(F) || return ω, it, ReturnCode.Success   # D exactly singular: ω is a root
        δ = inv(tr(F \ dD))
        isfinite(δ) || return oftype(ω, NaN), it, ReturnCode.Failure
        step = abs(δ) / max(abs(ω - δ), 1)
        (step <= atol || (step <= sqrt(atol) && step > prev / 2)) && return ω - δ, it, ReturnCode.Success
        it > 1 && (λ = min(1, abs(λ * δ₀) / 2abs(δ - (1 - λ) * δ₀)))
        ω -= λ * δ
        prev, δ₀ = step, δ
    end
    return oftype(ω, NaN), maxiter, ReturnCode.MaxIters
end

# A root of BOHH's own rational relation from ω0.
function seedsolver(prob, alg::BOHH; atol = 1.0e-12, maxiter = 50)
    T = realtype(prob)
    params = prepare(alg, prob.plasma.populations, T)
    return function (k, ω0)
        ω, nevals, retcode = @no_escape begin
            M = _operator(_temp_array, params, k, alg, Complex{T})
            r = size(M, 1) - 6
            E, Bf = (r + 1):(r + 3), (r + 4):(r + 6)
            _newton_det(M, r, M[E, Bf] * M[Bf, E], complex(ω0), atol, maxiter)
        end
        return (; omega = ω, resid = nothing, nevals, retcode)
    end
end
