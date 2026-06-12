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
    wc::T                       # Cyclotron frequency
    wp::T                      # Plasma frequency
    ρc::T                     # Cyclotron radius
    vtz::T                      # Parallel thermal velocity
    vtp::T                      # Perpendicular thermal velocity
    vdz::T                     # Parallel drift velocity
    vdr::T
    a0lm::Matrix{T}             # Normalized Hermite-Hermite coefficients
end

_convert_hh_param(param::HHSolverParam{T}, ::Type{T}) where {T} = param
function _convert_hh_param(param::HHSolverParam, ::Type{T}) where {T}
    return HHSolverParam{T}(
        T(param.wc), T(param.wp), T(param.ρc), T(param.vtz), T(param.vtp),
        T(param.vdz), T(param.vdr), Matrix{T}(param.a0lm),
    )
end

# Susceptibility is linear in a0lm.
# Negative density is realized as |n| with negated coefficients.
_a0lm(s::Maxwellian) = sign(s.n) * sqrt(eltype(s)(π)) * ones(eltype(s), 1, 1)
function _a0lm(vdf::BiKappa2)
    data = gen_fv2d(vdf)
    return hermite_expansion(data.fv, data.vz, data.vx, data.vtz, data.vtx)
end

_vtz(s) = sqrt(2 * kb * temperature(s.Tz) / s.m)
_vtp(s) = sqrt(2 * kb * temperature(s.Tp) / s.m)

function HHSolverParam(species, B0; a0lm = _a0lm(species))
    T = float(promote_type(typeof(B0), eltype(a0lm)))
    # Compute derived quantities for each species
    q = T(species.q)
    m = T(species.m)
    vtzs = T(_vtz(species))
    vtp = T(_vtp(species))
    wp = plasma_frequency(q, T(abs(species.n)), m)
    wc = T(B0) * q / m
    ρc = vtp / sqrt(T(2)) / wc
    vdz = T(species.vdz) * T(c0)
    vdr = T(species.vdr) * T(c0)
    return HHSolverParam{T}(wc, wp, ρc, vtzs, vtp, vdz, vdr, Matrix{T}(a0lm))
end

# Direct constructor from SI quantities.
# NB: `vtz`, `vtp`, `vdz`, `vdr` are in m/s here
# unlike distribution structs (`Maxwellian`, `BiKappa`) which use units of c.
function HHSolverParam(q, m, n, B0, vtz, vtp, vdz, vdr, a0lm)
    T = float(promote_type(map(typeof, (q, m, n, B0, vtz, vtp, vdz, vdr))..., eltype(a0lm)))
    wc = T(B0) * T(q) / T(m)
    wp = plasma_frequency(T(q), T(abs(n)), T(m))
    ρc = T(vtp) / sqrt(T(2)) / wc
    return HHSolverParam{T}(
        wc, wp, ρc, T(vtz), T(vtp), T(vdz), T(vdr), sign(n) * Matrix{T}(a0lm),
    )
end

# Precompute a0lm-weighted m-sums (independent of j) used inside the (n, j)-loop.
# Γₙ is the combined integral buffer from `_compute_integral_matrices!` with the
# component layout (A,p=1; A,p=2; B,p=1; B,p=2; C,p=1; C,p=2) along rows. For each
# z-Hermite order l and component k ∈ {A=1, B=2, C=3}:
#   Tₗ[k,l+1] = Σᵣ a0lm[l+1,r+1] Σₘ C[r+1,m+1] · (2Γₙ[2k-1,m+2] - mΓₙ[2k-1,m])
#   Uₗ[k,l+1] = Σᵣ a0lm[l+1,r+1] Σₘ C[r+1,m+1] · Γₙ[2k,m+1]
function _precompute_lm_sums!(Tₗ, Uₗ, a0lm, Cx, Γₙ)
    l_max, r_max = size(a0lm) .- 1
    @inbounds for l in 0:l_max
        ta = ua = tb = ub = tc = uc = zero(eltype(Tₗ))
        for r in 0:r_max
            a = a0lm[l + 1, r + 1]
            for m in 0:r
                c = Cx[r + 1, m + 1]
                dA = 2 * Γₙ[1, m + 2] - (m == 0 ? zero(eltype(Γₙ)) : m * Γₙ[1, m])
                dB = 2 * Γₙ[3, m + 2] - (m == 0 ? zero(eltype(Γₙ)) : m * Γₙ[3, m])
                dC = 2 * Γₙ[5, m + 2] - (m == 0 ? zero(eltype(Γₙ)) : m * Γₙ[5, m])
                ac = a * c
                ta += ac * dA
                ua += ac * Γₙ[2, m + 1]
                tb += ac * dB
                ub += ac * Γₙ[4, m + 1]
                tc += ac * dC
                uc += ac * Γₙ[6, m + 1]
            end
        end
        Tₗ[1, l + 1] = ta; Tₗ[2, l + 1] = tb; Tₗ[3, l + 1] = tc
        Uₗ[1, l + 1] = ua; Uₗ[2, l + 1] = ub; Uₗ[3, l + 1] = uc
    end
    return nothing
end

function _precompute_z_integrals!(Iₗ, DIₗ, EIₗ, Cz, Ils)
    l_max = length(Iₗ) - 1
    @inbounds for l in 0:l_max
        i = di = ei = zero(eltype(Iₗ))
        for k in 0:l
            c = Cz[l + 1, k + 1]
            i += c * Ils[k + 1]
            di += c * (2 * Ils[k + 2] - (k == 0 ? zero(eltype(Ils)) : k * Ils[k]))
            ei += c * (2 * Ils[k + 3] - k * Ils[k + 1])
        end
        Iₗ[l + 1] = i
        DIₗ[l + 1] = di
        EIₗ[l + 1] = ei
    end
    return nothing
end

function _precompute_z_poles!(Zₗ, DZₗ, czj)
    CT = eltype(Zₗ)
    T = _realtype(CT)
    l_max = size(Zₗ, 2) - 1
    p0 = inv(sqrt(sqrt(T(π))))
    @inbounds for j in axes(Zₗ, 1)
        c = czj[j]
        Zₗ[j, 1] = p0
        if l_max >= 1
            Zₗ[j, 2] = 2 * c * p0
            for l in 1:(l_max - 1)
                Zₗ[j, l + 2] = (2 * c / sqrt(T(l + 1))) * Zₗ[j, l + 1] -
                               sqrt(T(l) / T(l + 1)) * Zₗ[j, l]
            end
        end
        DZₗ[j, 1] = 2 * c * Zₗ[j, 1]
        for l in 1:l_max
            DZₗ[j, l + 1] = 2 * c * Zₗ[j, l + 1] - 2 * sqrt(T(l)) * Zₗ[j, l]
        end
    end
    return nothing
end

function _assemble_species_core!(
        M, Γₙ, Ils, Iₗ, DIₗ, EIₗ, Tₗ, Uₗ, Cz, Cx,
        snj, kx, kz, SNJ1, SNJ3,
        as, cj, bzj, Zₗ, DZₗ,
        wc, wp2, vtz, vtp, vdz, d, R, a0lm,
        N, J,
    )
    l_max, m_max = size(a0lm) .- 1
    z = zero(eltype(M))
    b11 = b12 = b13 = b22 = b23 = b33 = z

    vr = vtp / vtz
    dr = vdz / vtz
    vr2 = vr * vr
    kzvtz = kz * vtz
    kzvtp = kz * vtp
    vtzp = vtz / vtp

    Ils .= funIn.(Ref(eltype(Ils)), 0:(l_max + 2))
    _precompute_z_integrals!(Iₗ, DIₗ, EIₗ, Cz, Ils)

    for n in -N:N
        Γₙ .= 0
        nw_c = n * wc
        nwkp = nw_c / (kx * vtp)
        _compute_integral_matrices!(Γₙ, n, as, d, m_max)
        Γₙ .*= 2 / R

        _precompute_lm_sums!(Tₗ, Uₗ, a0lm, Cx, Γₙ)

        # j == 1 contributes per-n terms to species coefficients b11, b22, b33.
        # Lifted out of the j-loop since none of its factors depend on j.
        s11_3 = s22_3 = s33_3 = zero(eltype(Tₗ))
        @inbounds for l in 0:l_max
            s11_3 += Iₗ[l + 1] * Tₗ[1, l + 1]
            s22_3 += Iₗ[l + 1] * Tₗ[3, l + 1]
            s33_3 += (dr * DIₗ[l + 1] + EIₗ[l + 1]) * Uₗ[1, l + 1]
        end
        b11 -= wp2 * nwkp^2 * s11_3
        b22 -= wp2 * s22_3
        b33 -= wp2 * s33_3

        for j in 1:J
            snj += 1
            cⱼ = cj[j]
            cnj = cⱼ * kzvtz + kz * vdz + nw_c

            s11_1 = s11_2 = z
            s12_1 = s12_2 = z
            s22_1 = s22_2 = z
            s13_1 = s13_2 = z
            s23_1 = s23_2 = z
            s33_1 = s33_2 = z
            @inbounds for l in 0:l_max
                zₗ = Zₗ[j, l + 1]
                dzₗ = DZₗ[j, l + 1]
                cd = cⱼ + dr
                cd2 = cd * cd

                tA = Tₗ[1, l + 1]; uA = Uₗ[1, l + 1]
                tB = Tₗ[2, l + 1]; uB = Uₗ[2, l + 1]
                tC = Tₗ[3, l + 1]; uC = Uₗ[3, l + 1]

                s11_1 += zₗ * tA;      s11_2 += dzₗ * uA
                s12_1 += zₗ * tB;      s12_2 += dzₗ * uB
                s22_1 += zₗ * tC;      s22_2 += dzₗ * uC
                s13_1 += cd * zₗ * tA; s13_2 += cd * dzₗ * uA
                s23_1 += cd * zₗ * tB; s23_2 += cd * dzₗ * uB
                s33_1 += cd2 * zₗ * tA; s33_2 += cd2 * dzₗ * uA
            end

            tmp = wp2 * bzj[j] / cnj
            p11 = nwkp^2 * (nw_c * s11_1 + kzvtz * vr2 * s11_2) * tmp
            p12 = 1im * nwkp * (nw_c * s12_1 + kzvtz * vr2 * s12_2) * tmp
            p22 = (nw_c * s22_1 + kzvtz * vr2 * s22_2) * tmp
            p13 = nwkp * (vtzp * nw_c * s13_1 + kzvtp * s13_2) * tmp
            p23 = -1im * (vtzp * nw_c * s23_1 + kzvtp * s23_2) * tmp
            p33 = vtzp * (vtzp * nw_c * s33_1 + kzvtp * s33_2) * tmp

            jjx = snj
            jjy = snj + SNJ1
            jjz = snj + 2 * SNJ1

            M[jjx, jjx] = cnj
            M[jjx, (SNJ3 + 1):(SNJ3 + 3)] .= (p11, p12, p13)
            M[jjy, jjy] = cnj
            M[jjy, (SNJ3 + 1):(SNJ3 + 3)] .= (-p12, p22, p23)   # p21 = -p12
            M[jjz, jjz] = cnj
            M[jjz, (SNJ3 + 1):(SNJ3 + 3)] .= (p13, -p23, p33)   # p31 = p13, p32 = -p23

            b11 -= p11
            b12 -= p12
            b13 -= p13
            b22 -= p22
            b23 -= p23
            b33 -= p33
        end
    end
    return snj, (b11, b12, b13, b22, b23, b33)
end

function _assemble_species!(M, param::HHSolverParam{T}, snj, kx, kz, SNJ1, SNJ3, czj, bzj, Zₗ, DZₗ, N, J) where {T}
    as = kx * param.ρc * sqrt(T(2)) # Perpendicular wavenumber parameter
    vtp = param.vtp
    d = param.vdr / vtp
    R = exp(-d^2) + sqrt(T(π)) * d * erfc(-d) # Normalization parameters
    wp2 = param.wp^2
    a0lm = param.a0lm
    l_max, m_max = size(a0lm) .- 1
    n_max = max(l_max, m_max)

    return @no_escape begin
        Γₙ = @temp_array(T, 6, m_max + 3)
        Ils = @temp_array(T, l_max + 3)
        Iₗ = @temp_array(T, l_max + 1)
        DIₗ = @temp_array(T, l_max + 1)
        EIₗ = @temp_array(T, l_max + 1)
        Tₗ = @temp_array(T, 3, l_max + 1)
        Uₗ = @temp_array(T, 3, l_max + 1)
        C = @temp_array(T, n_max + 1, n_max + 1)
        _fill_hermite_coefficients!(C)
        Cz = @view C[1:(l_max + 1), 1:(l_max + 1)]
        Cx = @view C[1:(m_max + 1), 1:(m_max + 1)]
        _assemble_species_core!(
            M, Γₙ, Ils, Iₗ, DIₗ, EIₗ, Tₗ, Uₗ, Cz, Cx,
            snj, kx, kz, SNJ1, SNJ3, as, czj, bzj, Zₗ, DZₗ,
            param.wc, wp2, param.vtz, vtp, param.vdz, d, R, a0lm,
            N, J
        )
    end
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

prepare(::BOHH, species, B0, ::Type{T}) where {T} = map(species) do sp
    sp isa HHSolverParam ? _convert_hh_param(sp, T) : _convert_hh_param(HHSolverParam(sp, B0), T)
end

function dispersion_matrix!(M, pb::DispersionProblem, alg::BOHH; c2 = c0^2)
    T = _realtype(eltype(M))
    kx, kz = pb.kx, pb.kz
    N = alg.N
    params = prepare(alg, pb.species, pb.B0, T)
    S = length(params)
    (; J, bzj, czj) = get_jpole_coefficients(alg.J, T)

    # Handle singularities
    kx = iszero(kx) ? T(1.0e-30) : T(kx)
    kz = T(kz)
    # Compute matrix dimensions
    SNJ = S * (2 * N + 1) * J
    SNJ1 = SNJ + S
    SNJ3 = 3 * SNJ1

    # Adjust czj for kz sign
    kz < 0 && (czj = -czj)

    max_lmax = maximum(p -> size(p.a0lm, 1) - 1, params)

    @no_escape begin
        Zₗ = @temp_array(Complex{T}, alg.J, max_lmax + 1)
        DZₗ = @temp_array(Complex{T}, alg.J, max_lmax + 1)
        _precompute_z_poles!(Zₗ, DZₗ, czj)

        snj = 0
        for s in 1:S
            param = params[s]
            snj, (b11, b12, b13, b22, b23, b33) = _assemble_species!(
                M, param,
                snj, kx, kz,
                SNJ1, SNJ3,
                czj, bzj, Zₗ, DZₗ,
                N, J,
            )

            M[SNJ + s, SNJ3 + 1] = b11
            M[SNJ + s, SNJ3 + 2] = b12
            M[SNJ + s, SNJ3 + 3] = b13
            M[SNJ + SNJ1 + s, SNJ3 + 1] = -b12 # b21
            M[SNJ + SNJ1 + s, SNJ3 + 2] = b22
            M[SNJ + SNJ1 + s, SNJ3 + 3] = b23
            M[SNJ + 2SNJ1 + s, SNJ3 + 1] = b13 # b31
            M[SNJ + 2SNJ1 + s, SNJ3 + 2] = -b23 # b32
            M[SNJ + 2SNJ1 + s, SNJ3 + 3] = b33
        end
    end

    # E(J) coupling: J_xyz = j_xyz + sum(v_snj_xyz)
    M[SNJ3 + 1, 1:SNJ1] .= -1.0
    M[SNJ3 + 2, (SNJ1 + 1):2SNJ1] .= -1.0
    M[SNJ3 + 3, (2SNJ1 + 1):3SNJ1] .= -1.0
    _B_E_part!(M, SNJ3, kx, kz; c2)
    return M
end

function _B_E_part!(M, idx, kx, kz; c2 = c0^2)
    # E(B) coupling: Maxwell's equations
    M[idx + 1, idx + 5] = c2 * kz
    M[idx + 2, idx + 4] = -c2 * kz
    M[idx + 2, idx + 6] = c2 * kx
    M[idx + 3, idx + 5] = -c2 * kx

    # B(E) coupling: Faraday's law
    M[idx + 4, idx + 2] = -kz
    M[idx + 5, idx + 1] = kz
    M[idx + 5, idx + 3] = -kx
    M[idx + 6, idx + 2] = kx
    return M
end
