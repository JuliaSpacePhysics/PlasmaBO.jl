# Matrix eigenvalue solver for plasma dispersion relations
struct HHSolverParam{T}
    wc::T                       # Cyclotron frequency
    wp::T                      # Plasma frequency
    ρc::T                     # Cyclotron radius
    vtz::T                      # Parallel thermal velocity
    vtp::T                      # Perpendicular thermal velocity
    vdz::T                     # Parallel drift velocity
    vdr::T
    aslm::Matrix{T}             # Hermite expansion coefficients
end

HHSolverParam(param::HHSolverParam, B0) = param

# The coefficients a_{s,lm}
_alm(::Maxwellian{T}) where {T} = ones(T, 1, 1)

function HHSolverParam(species, B0)
    T = Float64
    # Compute derived quantities for each species
    q = charge(species)
    m = mass(species)
    Tz = temperature(species.Tz)  # eV -> K
    Tp = temperature(species.Tp)  # eV -> K

    vtzs = sqrt(2 * kb * Tz / m)
    vtp = sqrt(2 * kb * Tp / m)

    wp = plasma_frequency(q, species.n, m)
    wc = B0 * q / m
    ρc = sqrt(kb * Tp / m) / wc
    vdz = species.vdz * c0
    vdr = species.vdr * c0
    return HHSolverParam{T}(wc, wp, ρc, vtzs, vtp, vdz, vdr, _alm(species))
end

function HHSolverParam(q, m, n, B0, vtz, vtp, vdz, vdr, alm)
    T = Float64
    wc = B0 * q / m
    wp = plasma_frequency(q, n, m)
    ρc = vtp / sqrt(2) / wc
    return HHSolverParam{T}(wc, wp, ρc, vtz, vtp, vdz, vdr, alm)
end

function _assemble_species!(
        M,
        snj,
        kx, kz, SNJ1, SNJ3,
        as,
        cj, bzj, cj_l,
        wc, wp2, vtz, vtp, vdz, d, R, alm,
        N, J,
    )
    l_max, m_max = size(alm) .- 1

    b11, b12, b13, b21, b22, b23, b31, b32, b33 = ntuple(_ -> 0.0im, 9)
    Aₙ = zeros(m_max + 4, 2) # Γ_{a,n,m,p}
    Bₙ = zeros(m_max + 4, 2) # Γ_{b,n,m,p}
    Cₙ = zeros(m_max + 4, 2) # Γ_{c,n,m,p}

    vr = vtp / vtz
    dr = vdz / vtz

    Ils = funIn.(0:(l_max + 2))
    for n in -N:N
        Aₙ .= 0
        Bₙ .= 0
        Cₙ .= 0
        nw_c = n * wc
        nwkp = nw_c / (kx * vtp)

        for m in 0:(m_max + 2)
            idm = m + 2
            Aₙ[idm, 1] = 2.0 / R * funAn(n, as, d, m, 0)
            Aₙ[idm, 2] = 2.0 / R * funAn(n, as, d, m, 1)
            Bₙ[idm, 1] = 2.0 / R * funBn(n, as, d, m, 1)
            Bₙ[idm, 2] = 2.0 / R * funBn(n, as, d, m, 2)
            Cₙ[idm, 1] = 2.0 / R * funCn(n, as, d, m, 2)
            Cₙ[idm, 2] = 2.0 / R * funCn(n, as, d, m, 3)
        end

        for j in 1:J
            snj += 1
            cⱼ = cj[j]

            cnj = cⱼ * kz * vtz + kz * vdz + nw_c

            sum11tmp1, sum11tmp2, sum11tmp3 = 0.0im, 0.0im, 0.0im
            sum12tmp1, sum12tmp2 = 0.0im, 0.0im
            sum22tmp1, sum22tmp2, sum22tmp3 = 0.0im, 0.0im, 0.0im
            sum13tmp1, sum13tmp2 = 0.0im, 0.0im
            sum23tmp1, sum23tmp2 = 0.0im, 0.0im
            sum32tmp1, sum32tmp2 = 0.0im, 0.0im
            sum33tmp1, sum33tmp2, sum33tmp3 = 0.0im, 0.0im, 0.0im

            # @tullio sum11tmp1 := aslm[l + 1, m + 1] * czjj_l[$j, l + 2] * (2 * Ans[m + 1 + 2, 1] - m * Ans[m - 1 + 2, 1]) threads = false avx = false tensor = false
            # @tullio sum11tmp2 := aslm[l + 1, m + 1] * (2 * czjj_l[$j, l + 1 + 2] - l * czjj_l[$j, l + 2]) * Ans[m + 2, 2] threads = false avx = false tensor = false

            for l in 0:l_max
                cˡ = cj_l[j, l + 2]
                cˡ⁺¹ = cⱼ * cˡ
                cˡ⁺² = cⱼ * cˡ⁺¹
                cˡ⁺³ = cⱼ * cˡ⁺²
                cˡ⁻¹ = l >= 1 ? cj_l[j, l - 1 + 2] : 0.0im
                dZl = 2 * cˡ⁺¹ - l * cˡ⁻¹
                dZlc = dZl * cⱼ

                for m in 0:m_max
                    aₗₘ = alm[l + 1, m + 1] # aₗₘ
                    idx_m = m + 2
                    idx_mp1 = idx_m + 1
                    idx_mm1 = idx_m - 1

                    dAm = 2 * Aₙ[idx_mp1, 1] - m * Aₙ[idx_mm1, 1]
                    dBm = 2 * Bₙ[idx_mp1, 1] - m * Bₙ[idx_mm1, 1]
                    dCm = 2 * Cₙ[idx_mp1, 1] - m * Cₙ[idx_mm1, 1]

                    sum11tmp1 += aₗₘ * cˡ * dAm
                    sum11tmp2 += aₗₘ * dZl * Aₙ[idx_m, 2]

                    sum12tmp1 += aₗₘ * cˡ * dBm
                    sum12tmp2 += aₗₘ * dZl * Bₙ[idx_m, 2]

                    sum22tmp1 += aₗₘ * cˡ * dCm
                    sum22tmp2 += aₗₘ * dZl * Cₙ[idx_m, 2]

                    sum13tmp1 += aₗₘ * (dr * cˡ + cˡ⁺¹) * dAm
                    sum13tmp2 += aₗₘ * (dZlc + dr * dZl) * Aₙ[idx_m, 2]

                    sum23tmp1 += aₗₘ * (dr * cˡ + cˡ⁺¹) * dBm
                    sum23tmp2 += aₗₘ * (dZlc + dr * dZl) * Bₙ[idx_m, 2]

                    sum32tmp1 += aₗₘ * (dr * cˡ + cˡ⁺¹) * dBm
                    sum32tmp2 += aₗₘ * (dZlc + dr * dZl) * Bₙ[idx_m, 2]

                    sum33tmp1 += aₗₘ * (dr^2 * cˡ + 2 * dr * cˡ⁺¹ + cˡ⁺²) * dAm
                    sum33tmp2 += aₗₘ * (dr^2 * dZl + 2 * dr * dZlc + (2 * cˡ⁺³ - l * cˡ⁺¹)) * Aₙ[idx_m, 2]

                    if j == 1
                        Iₗ, Iₗ₊₁, Iₗ₊₂ = Ils[l + 1], Ils[l + 2], Ils[l + 3]
                        Iₗ₋₁ = l >= 1 ? Ils[l] : 0.0
                        sum11tmp3 += aₗₘ * Iₗ * dAm
                        sum22tmp3 += aₗₘ * Iₗ * dCm
                        sum33tmp3 += aₗₘ * (dr * (2 * Iₗ₊₁ - l * Iₗ₋₁) + (2 * Iₗ₊₂ - l * Iₗ)) * Aₙ[idx_m, 2]
                    end
                end
            end

            if j == 1
                nwkp = nw_c / (kx * vtp)
                b11 -= wp2 * nwkp^2 * sum11tmp3
                b22 -= wp2 * sum22tmp3
                b33 -= wp2 * sum33tmp3
            end

            tmp = wp2 * bzj[j] / cnj

            kzvtz = kz * vtz
            kzvtp = kz * vtp
            vr2 = vr * vr

            p11snj = nwkp^2 * (nw_c * sum11tmp1 + kzvtz * vr2 * sum11tmp2) * tmp
            p12snj = 1im * nwkp * (nw_c * sum12tmp1 + kzvtz * vr2 * sum12tmp2) * tmp
            p21snj = -p12snj
            p22snj = (nw_c * sum22tmp1 + kzvtz * vr2 * sum22tmp2) * tmp
            p13snj = nwkp * ((vtz / vtp) * nw_c * sum13tmp1 + kzvtp * sum13tmp2) * tmp
            p31snj = p13snj
            p23snj = -1im * ((vtz / vtp) * nw_c * sum23tmp1 + kzvtp * sum23tmp2) * tmp
            p32snj = 1im * ((vtz / vtp) * nw_c * sum32tmp1 + kzvtp * sum32tmp2) * tmp
            p33snj = (vtz / vtp) * ((vtz / vtp) * nw_c * sum33tmp1 + kzvtp * sum33tmp2) * tmp

            jjx = snj + 0 * SNJ1
            jjy = snj + 1 * SNJ1
            jjz = snj + 2 * SNJ1

            # # v_snjx equation: ω v_snjx = c_snj v_snjx + b_snj11 E_x + b_snj12 E_y + b_snj13 E_z
            M[jjx, jjx] = cnj
            M[jjx, (SNJ3 + 1):(SNJ3 + 3)] .= (p11snj, p12snj, p13snj)

            # v_snjy equation
            M[jjy, jjy] = cnj
            M[jjy, (SNJ3 + 1):(SNJ3 + 3)] .= (p21snj, p22snj, p23snj)

            # v_snjz equation
            M[jjz, jjz] = cnj
            M[jjz, (SNJ3 + 1):(SNJ3 + 3)] .= (p31snj, p32snj, p33snj)

            b11 -= p11snj
            b12 -= p12snj
            b21 -= p21snj
            b22 -= p22snj
            b13 -= p13snj
            b31 -= p31snj
            b23 -= p23snj
            b32 -= p32snj
            b33 -= p33snj
        end
    end
    return snj, (b11, b12, b13, b21, b22, b23, b31, b32, b33)
end

function _assemble_species!(
        M, param,
        snj, kx, kz,
        SNJ1, SNJ3,
        czj, bzj,
        czj_l, N,
        J,
    )
    as = kx * param.ρc * sqrt(2) # Perpendicular wavenumber parameter
    vtp = param.vtp
    d = param.vdr / vtp
    R = exp(-d^2) + sqrt(π) * d * erfc(-d) # Normalization parameters
    wp2 = param.wp^2
    return _assemble_species!(
        M, snj, kx, kz, SNJ1, SNJ3, as, czj, bzj, czj_l,
        param.wc, wp2, param.vtz, vtp, param.vdz, d, R, param.aslm,
        N, J
    )
end

"""
    solve_dispersion_matrix(params, kx, kz; J=8)

Solve the kinetic dispersion relation using the matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

This method transforms the dispersion relation into a matrix eigenvalue problem
using J-pole approximation for the plasma dispersion function, allowing
simultaneous computation of all wave modes.

- `kx`: Perpendicular wave vector component (m⁻¹)
- `kz`: Parallel wave vector component (m⁻¹)
- `J`: Number of poles for Z-function approximation (default: 8)

See also: [`get_jpole_coefficients`](@ref)
"""
function solve_dispersion_matrix end

# State Vector Organization
# The state vector has this structure:
# Indices 1 to SNJ1:           v_snj_x components + species j_x
# Indices SNJ1+1 to 2*SNJ1:    v_snj_y components + species j_y
# Indices 2*SNJ1+1 to 3*SNJ1:  v_snj_z components + species j_z
# Indices SNJ3+1 to SNJ3+6:    E_x, E_y, E_z, B_x, B_y, B_z
# where SNJ1 = SNJ + S (SNJ pole velocities + S species auxiliary j's)

# Compute matrix dimensions
_size(S::Int, N, J) = 3 * (S * (2 * N + 1) * J + S) + 6
_size(species, N, J) = _size(length(species), N, J)

function build_dispersion_matrix(species, args...; N = 2, J = 8, kw...)
    NN = _size(species, N, J)
    M = zeros(ComplexF64, NN, NN)
    return build_dispersion_matrix!(M, species, args...; N, J, kw...)
end

function build_dispersion_matrix!(M, params, kx, kz; N = 2, J = 8, c2 = c0^2)
    S = length(params)
    (; J, bzj, czj) = get_jpole_coefficients(J)

    # Handle singularities
    kx = kx == 0.0 ? 1.0e-30 : kx
    # Compute matrix dimensions
    SNJ = S * (2 * N + 1) * J
    SNJ1 = SNJ + S
    SNJ3 = 3 * SNJ1

    # Adjust czj for kz sign
    kz < 0 && (czj = -czj)

    # Precompute czj^l for all l values
    lsmax = [size(params[s].aslm, 1) - 1 for s in 1:S]
    max_lsmax = maximum(lsmax)
    czj_l = zeros(ComplexF64, J, max_lsmax + 5)
    for l in 0:(max_lsmax + 3)
        czj_l[:, l + 2] .= czj .^ l
    end

    snj = 0
    for s in 1:S
        param = params[s]
        snj, (b11, b12, b13, b21, b22, b23, b31, b32, b33) = _assemble_species!(
            M, param,
            snj, kx, kz,
            SNJ1, SNJ3,
            czj, bzj, czj_l,
            N, J,
        )

        M[SNJ + s, SNJ3 + 1] = b11
        M[SNJ + s, SNJ3 + 2] = b12
        M[SNJ + s, SNJ3 + 3] = b13
        M[SNJ + SNJ1 + s, SNJ3 + 1] = b21
        M[SNJ + SNJ1 + s, SNJ3 + 2] = b22
        M[SNJ + SNJ1 + s, SNJ3 + 3] = b23
        M[SNJ + 2SNJ1 + s, SNJ3 + 1] = b31
        M[SNJ + 2SNJ1 + s, SNJ3 + 2] = b32
        M[SNJ + 2SNJ1 + s, SNJ3 + 3] = b33
    end


    # E(J) coupling: J_xyz = j_xyz + sum(v_snj_xyz)
    M[SNJ3 + 1, 1:SNJ1] .= -1.0
    M[SNJ3 + 2, (SNJ1 + 1):2SNJ1] .= -1.0
    M[SNJ3 + 3, (2SNJ1 + 1):3SNJ1] .= -1.0

    # E(B) coupling: Maxwell's equations
    M[SNJ3 + 1, SNJ3 + 5] = c2 * kz
    M[SNJ3 + 2, SNJ3 + 4] = -c2 * kz
    M[SNJ3 + 2, SNJ3 + 6] = c2 * kx
    M[SNJ3 + 3, SNJ3 + 5] = -c2 * kx

    # B(E) coupling: Faraday's law
    M[SNJ3 + 4, SNJ3 + 2] = -kz
    M[SNJ3 + 5, SNJ3 + 1] = kz
    M[SNJ3 + 5, SNJ3 + 3] = -kx
    M[SNJ3 + 6, SNJ3 + 2] = kx
    return M
end

function solve_dispersion_matrix(params, kx, kz; kw...)
    M = build_dispersion_matrix(params, kx, kz; kw...)
    return solve_with_threads(6) do # 6 threads is faster than 8 threads
        eigvals!(M)
    end
end

function solve_dispersion_matrix!(M, params, kx, kz; kw...)
    fill!(M, zero(eltype(M)))
    build_dispersion_matrix!(M, params, kx, kz; kw...)
    return solve_with_threads(6) do
        eigvals!(M)
    end
end

function solve_kinetic_dispersion(species, B0, kx, kz; kw...)
    params = HHSolverParam.(species, B0)
    return solve_dispersion_matrix(params, kx, kz; kw...)
end

function solve_kinetic_dispersion(species, B0, ks::AbstractVector, θ; N = 2, J = 8, kw...)
    NN = _size(species, N, J)
    M = zeros(ComplexF64, NN, NN)
    params = HHSolverParam.(species, B0)
    ωs = map(ks) do k
        kx, kz = k .* sincos(θ)
        solve_dispersion_matrix!(M, params, kx, kz; N, J, kw...)
    end
    return DispersionSolution(ks, ωs)
end
