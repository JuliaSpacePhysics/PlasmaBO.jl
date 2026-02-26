using BenchmarkTools
using PlasmaBO
using PlasmaBO: q, me, mp, qe

const SUITE = BenchmarkGroup()

# ==============================================================================
# Umeda 2012 Case
# ==============================================================================
SUITE["Umeda2012"] = BenchmarkGroup()

let
    B0 = 96.24e-9  # [Tesla]
    T = 51 # [eV]
    ring_beam = Maxwellian(:e, 1.0e5, T; vdz = 0.1, vdr = 0.05)
    background = Maxwellian(:e, 9.0e5, T)

    species = [ring_beam, background]

    # Compute normalization
    wce = abs(B0 * q / me)
    lambdaD = Debye_length(species)
    kn = 1 / lambdaD

    # Wave vector: k*λD = 0.03, θ = 40°
    k = 0.03 / lambdaD
    θ = deg2rad(40)
    kx = k * sin(θ)
    kz = k * cos(θ)

    SUITE["Umeda2012"]["solve"] = @benchmarkable solve($species, $B0, $kx, $kz; N = 6, J = 12)

    k_ranges = (0.01:0.05:0.3) .* kn
    SUITE["Umeda2012"]["scan"] = @benchmarkable solve($species, $B0, $k_ranges, $θ; N = 6)

    sol = solve(species, B0, k_ranges, θ; N = 6)
    initial_points = [
        (0.1 * kn, 0.3im * wce),
        (0.1 * kn, 0.1im * wce),
        (0.2 * kn, 0.25im * wce),
    ]
    SUITE["Umeda2012"]["track"] = @benchmarkable track($sol, $initial_points)
end

# ==============================================================================
# Fluid vs Kinetic
# ==============================================================================
SUITE["FluidVsKinetic"] = BenchmarkGroup()

let
    # Cold plasma parameters
    B0 = 1.0e-8  # [Tesla]
    n = 1.0e6    # [m^-3]
    T = 0.01     # [eV] - very cold

    # Wave vector
    k = 1.0  # [m^-1]
    θ = deg2rad(45)
    kx, kz = k .* sincos(θ)

    # Kinetic solver with cold Maxwellian
    e_vdf = Maxwellian(:e, n, T)
    kinetic_species = [e_vdf, Maxwellian(:p, n, T)]

    SUITE["FluidVsKinetic"]["kinetic"] = @benchmarkable solve($kinetic_species, $B0, $kx, $kz; N = 2, J = 8)

    # Fluid solver
    fluid_species = [e_vdf, FluidSpecies(:p, n, T)]
    SUITE["FluidVsKinetic"]["fluid"] = @benchmarkable solve($fluid_species, $B0, $kx, $kz, BOFluid)
end

# ==============================================================================
# Mirror Mode
# ==============================================================================
SUITE["MirrorMode"] = BenchmarkGroup()

let
    B0 = 100.0e-9
    θ = deg2rad(71)

    n = 1.0e6
    Tpara = 24840.0
    Tperp = 49680.0

    electron = Maxwellian(:e, n, Tpara)
    proton = Maxwellian(n, Tpara, Tperp)

    species = (proton, electron)
    params = HHSolverParam.(species, B0)
    ρᵢ = params[1].ρc

    # Selecting a single k for benchmark, from the range (0.005:0.04:0.7) ./ ρᵢ
    # Choosing a value in the middle
    k_val = 0.35 / ρᵢ
    kx, kz = k_val .* sincos(θ)

    SUITE["MirrorMode"]["solve"] = @benchmarkable solve($species, $B0, $kx, $kz)
end
