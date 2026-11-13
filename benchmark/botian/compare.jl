using PlasmaBO, LinearAlgebra, Printf
using PlasmaPopulations: Particle

BLAS.set_num_threads(1)  # as BO-T; (k, θ) sweeps still use Julia threads

# BO-T's SI constants (bo_setup.m); its electron is 5.447e-4 proton masses.
const qe, mp, c = 1.60217662e-19, 1.6726219e-27, 2.99792458e8
const proton, electron = Particle(qe, mp), Particle(-qe, 5.447e-4 * mp)
vth(p, T_eV) = sqrt(2 * T_eV * qe / p.m)

results(file) = joinpath(@__DIR__, "results", file)

"BO-T's wavevectors, eigenvalues at each, and total eigensolve seconds."
function load(name)
    ks, roots, seconds = Wavevector[], Vector{ComplexF64}[], 0.0
    for l in eachline(results("$name.csv"))
        i, kx, kz, t, re, im = parse.(Float64, split(l, ','))
        if i == length(ks)
            push!(ks, Wavevector(kx, kz)); push!(roots, ComplexF64[]); seconds += t
        end
        push!(roots[end], complex(re, im))
    end
    return ks, roots, seconds
end

# Each returns the plasma, the solver at BO-T's (N, J), ωn and BO-T's plot window [ωn].
function umeda()
    v = vth(electron, 51.0)
    plasma = Plasma(
        Population(electron, GaussianRing(; vr_perp = 0.05c, w_perp = v, vth_para = v, vd_para = 0.1c); n = 1.0e5),
        Population(electron, BiMaxwellian(; vth = v); n = 9.0e5);
        B0 = 96.24e-9,
    )
    return plasma, BOHH(; N = 8, J = 8), qe * 96.24e-9 / electron.m, (-8, 8, -2, 0.5)
end

function firehose()
    lines = collect(eachline(results("firehose_alm.csv")))
    w_perp, w_para = parse.(Float64, match(r"w_perp=(\S+) w_para=(\S+)", lines[1]).captures)
    alm = stack(parse.(Float64, split(l, ',')) for l in lines[2:end]; dims = 1)
    n, B0 = 5.0e19, 0.1
    plasma = Plasma(
        Population(proton, HermiteVDF(; w_perp, w_para, alm); n),
        Population(electron, BiMaxwellian(; vth = vth(electron, 496.683)); n);
        B0,
    )
    return plasma, BOHH(; N = 2, J = 24), qe * B0 / mp, (-1, 1, -0.05, 0.1)
end

timed(f) = (f(); minimum(@elapsed(f()) for _ in 1:3))

"Every root at BO-T's 1-D scan: timings, both root sets in BO-T's plot window [ωn], and the worst
distance of a PlasmaBO root from BO-T's, |Δω| / (ωn + |ω|) as in BO-T."
function scan(name, setup)
    ks, ref, t_ref = load(name)
    plasma, alg, ωn, (xmin, xmax, ymin, ymax) = setup()
    sol(k) = solve(DispersionProblem(plasma, k), alg).omega
    t = timed(() -> foreach(sol, ks))
    inwindow(ω) = xmin <= real(ω) <= xmax && ymin <= imag(ω) <= ymax
    k = [abs(complex(k.k_perp, k.k_para)) / abs(complex(ks[1].k_perp, ks[1].k_para)) for k in ks]
    ours = [filter(inwindow, sol(k) ./ ωn) for k in ks]
    theirs = [filter(inwindow, r ./ ωn) for r in ref]
    err = maximum(maximum((minimum(abs, r .- ω) / (1 + abs(ω)) for ω in o); init = 0.0) for (o, r) in zip(ours, ref ./ ωn))
    return (; alg, t_ref, t, err, k, ours, theirs)
end

"The ring-beam instability over a (k, θ) grid: every root, or Newton continuation of its branch."
function map_scan()
    ks, ref, t_ref = load("umeda_map")
    plasma, alg, ωn, _ = umeda()
    kλ, θ = (1:30) ./ 100, deg2rad.(5:5:85)
    kn = abs(complex(ks[1].k_perp, ks[1].k_para)) / kλ[1]
    sweep = AngleSweep(; k = kn .* kλ, theta = θ)
    t_survey = timed(() -> solve(DispersionProblem(plasma, sweep), alg))

    # Seed from one eigensolve, continue along θ, then along k from each θ.
    at = Wavevector(; k = sweep.k[12], theta = θ[1])
    ω0 = argmax(imag, solve(DispersionProblem(plasma, at), alg).omega)
    follow(ω, at, k, theta) = vec(solve(DispersionProblem(plasma, Seed(ω, at), AngleSweep(; k, theta)), Continuation(alg)).omega)
    function branch()
        seeds = follow(ω0, at, sweep.k[12], θ)
        return stack(follow(seeds[j], Wavevector(; k = sweep.k[12], theta = θ[j]), sweep.k, θ[j]) for j in eachindex(θ))
    end
    t_branch = timed(branch)
    γ = imag.(branch()) ./ ωn
    γ_ref = maximum(maximum(imag, r) for r in ref) / ωn
    return (; alg, n = length(ks), t_ref, t_survey, t_branch, kλ, θ = rad2deg.(θ), γ, γ_ref)
end

ring, fire, grid = scan("umeda", umeda), scan("firehose", firehose), map_scan()

println("Julia threads: ", Threads.nthreads(), "\n")
println("| task | PlasmaBO solver | BO-T | PlasmaBO | speedup | PlasmaBO root vs BO-T |")
println("|---|---|---|---|---|---|")
for (label, s) in (("ring-beam, all roots, $(length(ring.k)) k", ring), ("firehose, all roots, $(length(fire.k)) k", fire))
    @printf("| %s | `%s` | %.1f s | %.1f s | **%.1f×** | ≤ %.0e |\n", label, s.alg, s.t_ref, s.t, s.t_ref / s.t, s.err)
end
@printf("| ring-beam (k, θ) map, all roots, %d k | `%s`, %d threads | %.0f s | %.1f s | **%.0f×** | |\n",
    grid.n, grid.alg, Threads.nthreads(), grid.t_ref, grid.t_survey, grid.t_ref / grid.t_survey)
@printf("| ring-beam (k, θ) map, unstable branch | `Continuation(%s)` | %.0f s | %.2f s | **%.0f×** | max γ/Ωe %.6f vs %.6f |\n",
    grid.alg, grid.t_ref, grid.t_branch, grid.t_ref / grid.t_branch, maximum(grid.γ), grid.γ_ref)

using CairoMakie
fig = Figure(size = (1400, 430))
let ax = Axis(fig[1, 1]; title = "(a) wall time [s]", xscale = log10, yreversed = true,
        yticks = (1:4, ["ring-beam\nall roots", "firehose\nall roots", "(k, θ) map\nall roots", "(k, θ) map\none branch"]))
    ts = [(ring.t_ref, ring.t), (fire.t_ref, fire.t), (grid.t_ref, grid.t_survey), (grid.t_ref, grid.t_branch)]
    for (i, (a, b)) in enumerate(ts)
        barplot!(ax, [i - 0.2], [a]; direction = :x, width = 0.38, color = :gray70, fillto = 0.1)
        barplot!(ax, [i + 0.2], [b]; direction = :x, width = 0.38, color = :dodgerblue, fillto = 0.1)
        text!(ax, b, i + 0.2; text = @sprintf(" %.0f×", a / b), align = (:left, :center), fontsize = 13)
    end
    xlims!(ax, 0.1, 400)
    axislegend(ax, [PolyElement(color = :gray70), PolyElement(color = :dodgerblue)], ["BO-T", "PlasmaBO"]; position = :rb)
end
let ax = Axis(fig[1, 2]; title = "(b) ring-beam roots, BO-T (gray) and PlasmaBO (blue)", xlabel = "kλD", ylabel = "Im ω/Ωe")
    for (s, color, size) in ((ring.theirs, :gray70, 7), (ring.ours, :dodgerblue, 3))
        x = reduce(vcat, fill(0.01k, length(r)) for (k, r) in zip(ring.k, s))
        scatter!(ax, x, imag.(reduce(vcat, s)); color, markersize = size)
    end
    ylims!(ax, -1, 0.4)
end
let ax = Axis(fig[1, 3]; title = "(c) γ/Ωe of the branch continued from θ = 5°", xlabel = "kλD", ylabel = "θ [°]")
    hm = heatmap!(ax, grid.kλ, grid.θ, grid.γ; colormap = :inferno, colorrange = (0, maximum(grid.γ)))
    Colorbar(fig[1, 4], hm)
end
save(joinpath(@__DIR__, "benchmark.png"), fig; px_per_unit = 1.5)
