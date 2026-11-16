using BenchmarkTools
using PlasmaBO
using PlasmaPopulations: gyrofrequency, thermal_speed, debye_length, kappa_widths
using PlasmaPopulations: C_SI as c0

const SUITE = BenchmarkGroup()

let g = SUITE["Umeda2012"]
    B0, T, e = 96.24e-9, 51, Electron
    vth = thermal_speed(e, T)
    plasma = Plasma(
        Population(e, GaussianRing(; vr_perp = 0.05c0, w_perp = vth, vth_para = vth, vd_para = 0.1c0); n = 1.0e5),
        Population(e, BiMaxwellian; T, n = 9.0e5);
        B0,
    )
    wce = abs(gyrofrequency(e, B0))
    kn = 1 / debye_length(e, 1.0e6, T)
    θ = deg2rad(40)

    g["solve"] = @benchmarkable solve($(DispersionProblem(plasma, Wavevector(; k = 0.03kn, theta = θ))), BOHH(; N = 6, J = 12))
    sweep = DispersionProblem(plasma, AngleSweep(; k = collect(0.01:0.05:0.3) .* kn, theta = θ))
    g["scan"] = @benchmarkable solve($sweep, BOHH(; N = 6))
    survey = solve(sweep, BOHH(; N = 6))
    g["link"] = @benchmarkable link($survey)
    g["track"] = @benchmarkable track($survey, $((2.32 + 0.33im) * wce), $(Wavevector(; k = 0.1kn, theta = θ)))
end

let g = SUITE["MirrorMode"]
    B0, n, p = 100.0e-9, 1.0e6, Proton
    plasma = Plasma(
        Population(p, BiMaxwellian; T_perp = 49680.0, T_para = 24840.0, n),
        Population(Electron, BiMaxwellian; T = 24840.0, n);
        B0,
    )
    ρi = thermal_speed(p, 49680.0) / sqrt(2) / gyrofrequency(p, B0)
    problem = DispersionProblem(plasma, Wavevector(; k = 0.35 / ρi, theta = deg2rad(71)))
    g["BOHH"] = @benchmarkable solve($problem, BOHH())
    g["BOFluid"] = @benchmarkable solve($problem, BOFluid())
end

let g = SUITE["Cattaert2007"]
    B0, e = 1.0e-6, Electron
    vdf = ProductBiKappa(e; T = 2555.0, kappa_perp = 200.0, kappa_para = 1.0)
    w_perp, w_para = kappa_widths(vdf)
    plasma = Plasma(Population(e, vdf; n = 2.43e6); B0)
    kn = abs(gyrofrequency(e, B0)) / w_perp
    θ = deg2rad(30)

    g["HermiteVDF"] = @benchmarkable HermiteVDF($vdf; w_perp = $w_perp, w_para = $w_para, N_para = 16, N_perp = 16)
    k = Wavevector(; k = 0.1kn, theta = θ)
    g["matrix BOHH"] = @benchmarkable dispersion_matrix($plasma, $k, BOHH())
    g["matrix BOPBK"] = @benchmarkable dispersion_matrix($plasma, $k, BOPBK())
    g["scan"] = @benchmarkable solve($(DispersionProblem(plasma, AngleSweep(; k = collect(range(1.0e-4, 0.3; length = 80) .* kn), theta = θ))), BOPBK())
end
