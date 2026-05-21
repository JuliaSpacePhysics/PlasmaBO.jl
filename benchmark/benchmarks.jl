using BenchmarkTools
using PlasmaBO
using PlasmaBO: q, me, mp, qe, interpolate_complex

const SUITE = BenchmarkGroup()

let
    # PCHIP path: needs length(z_prev) ≥ 3; size matches a typical k-scan run.
    x_prev = collect(range(0.0, 1.0; length = 16))
    z_prev = ComplexF64.(2.0 .* x_prev .+ 1.0, 0.1 .* x_prev .+ 0.5 .* sin.(x_prev))
    x_new = 1.05  # slight extrapolation, as in real tracking
    SUITE["Track"]["interpolate_complex_pchip"] =
        @benchmarkable interpolate_complex($z_prev, $x_prev, $x_new)
end

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
    SUITE["Umeda2012"]["track"] = @benchmarkable track.($sol, $initial_points)
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

    species = [Maxwellian(:e, n, T), Maxwellian(:p, n, T)]
    SUITE["FluidVsKinetic"]["kinetic"] = @benchmarkable solve($species, $B0, $kx, $kz; N = 2, J = 8)
    SUITE["FluidVsKinetic"]["fluid"] = @benchmarkable solve($species, $B0, $kx, $kz, BOFluid)
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

SUITE["HermiteExpansion"] = BenchmarkGroup()

let
    bk = BiKappa2(:e, 1.0e6, 1.0, 200.0, 2555.0; sigma = 0.0)
    g = gen_fv2d(bk)
    SUITE["HermiteExpansion"]["BiKappa2"] = @benchmarkable hermite_expansion(
        $(g.fv), $(g.vz[:, 1]), $(g.vx[1, :]), $(g.vtz), $(g.vtx); Nz = 16, Nx = 16
    )
end

# ==============================================================================
# PBK Solver
# ==============================================================================
SUITE["PBKSolver"] = BenchmarkGroup()

let
    B0 = 1.0e-6          # [Tesla]
    θ = deg2rad(30)

    n = 2.43e6           # [m^-3]
    T = 2555.0           # [eV]
    κz = 1.0
    κx = 200.0
    electron = BiKappa2(:e, n, κz, κx, T; sigma = 0.0)

    wce = abs(B0 * q / me)
    ρce = electron.vtp / wce
    kn = 1 / ρce  # so that k/kn = kρce
    kρ_scan = range(1.0e-4, 0.3; length = 80)
    ks = kρ_scan .* kn

    SUITE["PBKSolver"]["solve"] = @benchmarkable solve($electron, $B0, $ks, $θ, BOPBK)
end
