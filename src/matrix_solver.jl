# Matrix eigenvalue solver for plasma dispersion relations
"""
    HHSolverParams{T}

Parameters for the matrix eigenvalue solver.
"""
struct HHSolverParams{T} <: AbstractSolverParams
    S::Int                      # Number of species
    c2::T                       # Speed of light squared
    wcs::Vector{T}              # Cyclotron frequencies
    wps2::Vector{T}             # Plasma frequencies squared
    rhocs::Vector{T}            # Cyclotron radii
    vtzs::Vector{T}             # Parallel thermal velocities
    vtps::Vector{T}             # Perpendicular thermal velocities
    vdsz::Vector{T}             # Parallel drift velocities
    ds::Vector{T}               # Ring beam drift parameters
    Rs::Vector{T}               # Normalization parameters
    aslm::Vector{Matrix{T}}     # Hermite expansion coefficients
    msmax::Vector{Int}          # Max perpendicular Hermite index in HH expansion
    lsmax::Vector{Int}          # Max parallel Hermite index in HH expansion
end


"""
    create_solver_params(species::Vector{Species}, B0, kx, kz)

Create HHSolverParams from a list of Species.

# Arguments
- `species`: Vector of Species objects
- `B0`: Magnetic field strength (Tesla)
- `kx`: Perpendicular wave vector (m⁻¹)
- `kz`: Parallel wave vector (m⁻¹)

See also: [`Species`](@ref), [`fDrHH`](@ref)
"""
function create_solver_params(species, B0)
    S = length(species)
    T = Float64

    c2 = c0^2

    # Compute derived quantities for each species
    qs = charge.(species)
    ms = mass.(species)
    ns = [s.n for s in species]
    Tzs = [s.Tz * q / kb for s in species]  # eV -> K
    Tps = [s.Tp * q / kb for s in species]

    vtzs = @. sqrt(2 * kb * Tzs / ms)
    vtps = @. sqrt(2 * kb * Tps / ms)

    wps = plasma_frequency.(qs, ns, ms)
    wps2 = wps .^ 2

    wcs = @. B0 * qs / ms
    rhocs = @. sqrt(kb * Tps / ms) / wcs

    vdsz = [sp.vdz * c0 for sp in species]
    ds = [species[i].vdr * c0 / vtps[i] for i in 1:S]

    As = [exp(-ds[i]^2) + sqrt(π) * ds[i] * erfc(-ds[i]) for i in 1:S]

    aslm = [s.aslm for s in species]
    msmax = [size(s.aslm, 2) - 1 for s in species]
    lsmax = [size(s.aslm, 1) - 1 for s in species]

    return HHSolverParams{T}(
        S, T(c2), T.(wcs), T.(wps2), T.(rhocs),
        T.(vtzs), T.(vtps), T.(vdsz),
        T.(ds), T.(As), aslm, msmax, lsmax
    )
end


function _assemble_species!(
        M,
        s,
        snj,
        kx, kz,
        S, SNJ1, SNJ3,
        as_s,
        cj, cj_l,
        wcₛ, wp2ₛ, vtzₛ, vtpₛ, vdzₛ, dₛ,
        R, # normalization factor
        N,
        aslm,
        m_max, l_max,
        bzj,
        J,
    )

    b11, b12, b13, b21, b22, b23, b31, b32, b33 = ntuple(_ -> 0.0im, 9)
    Aₙ = zeros(m_max + 4, 2) # Γ_{a,n,m,p}
    Bₙ = zeros(m_max + 4, 2) # Γ_{b,n,m,p}
    Cₙ = zeros(m_max + 4, 2) # Γ_{c,n,m,p}

    vr = vtpₛ / vtzₛ
    dr = vdzₛ / vtzₛ

    Ils = funIn.(0:(l_max + 2))
    for n in -N:N
        Aₙ .= 0
        Bₙ .= 0
        Cₙ .= 0
        nw_c = n * wcₛ
        nwkp = nw_c / (kx * vtpₛ)

        for m in 0:(m_max + 2)
            idm = m + 2
            Aₙ[idm, 1] = 2.0 / R * funAn(n, as_s, dₛ, m, 0)
            Aₙ[idm, 2] = 2.0 / R * funAn(n, as_s, dₛ, m, 1)
            Bₙ[idm, 1] = 2.0 / R * funBn(n, as_s, dₛ, m, 1)
            Bₙ[idm, 2] = 2.0 / R * funBn(n, as_s, dₛ, m, 2)
            Cₙ[idm, 1] = 2.0 / R * funCn(n, as_s, dₛ, m, 2)
            Cₙ[idm, 2] = 2.0 / R * funCn(n, as_s, dₛ, m, 3)
        end

        for j in 1:J
            snj += 1
            cⱼ = cj[j]

            cnj = cⱼ * kz * vtzₛ + kz * vdzₛ + nw_c

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
                    aₗₘ = aslm[l + 1, m + 1] # aₗₘ
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
                nwkp = nw_c / (kx * vtpₛ)
                b11 -= wp2ₛ * nwkp^2 * sum11tmp3
                b22 -= wp2ₛ * sum22tmp3
                b33 -= wp2ₛ * sum33tmp3
            end

            tmp = wp2ₛ * bzj[j] / cnj

            kzvtz = kz * vtzₛ
            kzvtp = kz * vtpₛ
            vr2 = vr * vr

            p11snj = nwkp^2 * (nw_c * sum11tmp1 + kzvtz * vr2 * sum11tmp2) * tmp
            p12snj = 1im * nwkp * (nw_c * sum12tmp1 + kzvtz * vr2 * sum12tmp2) * tmp
            p21snj = -p12snj
            p22snj = (nw_c * sum22tmp1 + kzvtz * vr2 * sum22tmp2) * tmp
            p13snj = nwkp * ((vtzₛ / vtpₛ) * nw_c * sum13tmp1 + kzvtp * sum13tmp2) * tmp
            p31snj = p13snj
            p23snj = -1im * ((vtzₛ / vtpₛ) * nw_c * sum23tmp1 + kzvtp * sum23tmp2) * tmp
            p32snj = 1im * ((vtzₛ / vtpₛ) * nw_c * sum32tmp1 + kzvtp * sum32tmp2) * tmp
            p33snj = (vtzₛ / vtpₛ) * ((vtzₛ / vtpₛ) * nw_c * sum33tmp1 + kzvtp * sum33tmp2) * tmp

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

    SNJ = SNJ1 - S
    M[SNJ + s, SNJ3 + 1] = b11
    M[SNJ + s, SNJ3 + 2] = b12
    M[SNJ + s, SNJ3 + 3] = b13
    M[SNJ + SNJ1 + s, SNJ3 + 1] = b21
    M[SNJ + SNJ1 + s, SNJ3 + 2] = b22
    M[SNJ + SNJ1 + s, SNJ3 + 3] = b23
    M[SNJ + 2SNJ1 + s, SNJ3 + 1] = b31
    M[SNJ + 2SNJ1 + s, SNJ3 + 2] = b32
    M[SNJ + 2SNJ1 + s, SNJ3 + 3] = b33

    return snj
end

"""
    solve_dispersion_matrix(params, kx, kz; J=8)

Solve the kinetic dispersion relation using the matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

This method transforms the dispersion relation into a matrix eigenvalue problem
using J-pole approximation for the plasma dispersion function, allowing
simultaneous computation of all wave modes.

# Arguments
- `params`: HHSolverParams or MatrixSolverParams with plasma parameters
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

function build_dispersion_matrix(params::HHSolverParams, kx, kz; N = 2, J = 8)
    (; S, c2, wcs, wps2, rhocs, vtzs, vtps, vdsz, ds, Rs, aslm, msmax, lsmax) = params
    (; J, bzj, czj) = get_jpole_coefficients(J)

    # Handle singularities
    kx = kx == 0.0 ? 1.0e-30 : kx
    # Perpendicular wavenumber parameter
    as = kx .* rhocs .* sqrt(2)

    # Compute matrix dimensions
    Ns = 2 .* fill(N, S) .+ 1  # Number of harmonics per species
    SNJ = S * (2 * N + 1) * J
    SNJ1 = SNJ + S
    SNJ3 = 3 * SNJ1
    NN = SNJ3 + 6

    # Adjust czj for kz sign
    kz < 0 && (czj = -czj)

    # Precompute czj^l for all l values
    max_lsmax = maximum(lsmax)
    czj_l = zeros(ComplexF64, J, max_lsmax + 5)
    for l in 0:(max_lsmax + 3)
        czj_l[:, l + 2] .= czj .^ l
    end

    # Initialize coefficient arrays
    # M = Matrix(sparse(I_idx, J_idx, V_val, NN, NN))
    M = zeros(ComplexF64, NN, NN)

    snj = 0
    for s in 1:S
        snj = _assemble_species!(
            M,
            s,
            snj,
            kx,
            kz,
            S,
            SNJ1,
            SNJ3,
            as[s],
            czj,
            czj_l,
            wcs[s],
            wps2[s],
            vtzs[s],
            vtps[s],
            vdsz[s],
            ds[s],
            Rs[s],
            N,
            aslm[s],
            msmax[s],
            lsmax[s],
            bzj,
            J,
        )
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

function solve_kinetic_dispersion(species, B0, kx, kz; kw...)
    params = create_solver_params(species, B0)
    return solve_dispersion_matrix(params, kx, kz; kw...)
end


function solve_kinetic_dispersion(species, B0, ks::AbstractVector, θ; kw...)
    ωs = map(ks) do k
        kx = k * sin(θ)
        kz = k * cos(θ)
        solve_kinetic_dispersion(species, B0, kx, kz; kw...)
    end
    return DispersionSolution(ks, ωs)
end
