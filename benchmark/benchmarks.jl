using BenchmarkTools
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed, debye_length
using PlasmaPopulations: E_SI as q, ME_SI as me, C_SI as c0, kappa_widths
using PlasmaBO: hermite_expansion, gen_fv2d

const SUITE = BenchmarkGroup()

cattaert_electron() = ProductBiKappa(Electron(); T = 2555.0, kappa_perp = 200.0, kappa_para = 1.0)

SUITE["Umeda2012"] = BenchmarkGroup()

let
    B0 = 96.24e-9  # [Tesla]
    T = 51 # [eV]
    e = Electron()
    plasma = Plasma(
        Population(e, (vth = thermal_speed(e, T); GaussianRing(; vr_perp = 0.05c0, w_perp = vth, vth_para = vth, vd_para = 0.1c0)); n = 1.0e5),
        Population(e, BiMaxwellian; T, n = 9.0e5);
        B0,
    )

    # Compute normalization
    wce = abs(B0 * q / me)
    lambdaD = debye_length(e, 1.0e6, T)
    kn = 1 / lambdaD

    # Wave vector: k*λD = 0.03, θ = 40°
    k = 0.03 / lambdaD
    θ = deg2rad(40)
    kx = k * sin(θ)
    kz = k * cos(θ)

    SUITE["Umeda2012"]["solve"] = @benchmarkable solve($(DispersionProblem(plasma, Wavevector(kx, kz))), BOHH(; N = 6, J = 12))

    k_ranges = (0.01:0.05:0.3) .* kn
    sweep = AngleSweep(; k = collect(k_ranges), theta = θ)
    SUITE["Umeda2012"]["scan"] = @benchmarkable solve($(DispersionProblem(plasma, sweep)), BOHH(; N = 6))
    ω0, k0 = (2.32 + 0.33im) * wce, Wavevector((0.1kn) .* sincos(θ)...)
    survey = solve(DispersionProblem(plasma, sweep), BOHH(; N = 6))
    SUITE["Umeda2012"]["link"] = @benchmarkable link($survey)
    SUITE["Umeda2012"]["track"] = @benchmarkable track($survey, $ω0, $k0)
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

    plasma = Plasma(Population(Electron(), BiMaxwellian; T, n), Population(Proton(), BiMaxwellian; T, n); B0)
    SUITE["FluidVsKinetic"]["kinetic"] = @benchmarkable solve($(DispersionProblem(plasma, Wavevector(kx, kz))), BOHH(; N = 2, J = 8))
    SUITE["FluidVsKinetic"]["fluid"] = @benchmarkable solve($(DispersionProblem(plasma, Wavevector(kx, kz))), BOFluid())
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

    p = Proton()
    plasma = Plasma(
        Population(p, BiMaxwellian; T_perp = Tperp, T_para = Tpara, n),
        Population(Electron(), BiMaxwellian; T = Tpara, n);
        B0,
    )
    ρᵢ = thermal_speed(p, Tperp) / sqrt(2) / gyrofrequency(p, B0)

    # Selecting a single k for benchmark, from the range (0.005:0.04:0.7) ./ ρᵢ
    # Choosing a value in the middle
    k_val = 0.35 / ρᵢ
    kx, kz = k_val .* sincos(θ)

    SUITE["MirrorMode"]["solve"] = @benchmarkable solve($(DispersionProblem(plasma, Wavevector(kx, kz))), BOHH())
end

SUITE["HermiteExpansion"] = BenchmarkGroup()

let
    g = gen_fv2d(cattaert_electron())
    SUITE["HermiteExpansion"]["ProductBiKappa"] = @benchmarkable hermite_expansion(
        $(g.fv), $(g.v_para), $(g.v_perp), $(g.w_para), $(g.w_perp); N_para = 16, N_perp = 16
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
    electron = Plasma(Population(Electron(), cattaert_electron(); n); B0)

    wce = abs(B0 * q / me)
    ρce = kappa_widths(cattaert_electron())[1] / wce
    kn = 1 / ρce  # so that k/kn = kρce
    kρ_scan = range(1.0e-4, 0.3; length = 80)
    ks = kρ_scan .* kn

    SUITE["PBKSolver"]["solve"] = @benchmarkable solve($(DispersionProblem(electron, AngleSweep(; k = collect(ks), theta = θ))), BOPBK())
end


# Matrix Construction - isolating matrix assembly from eigensolve that dominates `solve`.
SUITE["BuildMatrix"] = BenchmarkGroup()

let
    B0 = 1.0e-6
    electron = Plasma(Population(Electron(), cattaert_electron(); n = 2.43e6); B0)
    wce = abs(B0 * q / me)
    ρce = kappa_widths(cattaert_electron())[1] / wce
    kn = 1 / ρce
    θ = deg2rad(30)
    kx = 0.1 * kn * cos(θ)
    kz = 0.1 * kn * sin(θ)

    SUITE["BuildMatrix"]["HH"] = @benchmarkable dispersion_matrix($electron, $(Wavevector(kx, kz)), BOHH())
    SUITE["BuildMatrix"]["PBK"] = @benchmarkable dispersion_matrix($electron, $(Wavevector(kx, kz)), BOPBK())
end
