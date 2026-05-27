module PrecisionSolutionPlot

using CairoMakie
using DoubleFloats
using GenericLinearAlgebra
using LinearAlgebra
using PlasmaBO
using Printf
using SpecialFunctions

const B0 = 0.1
const THETA = deg2rad(45)
const DENSITY = 5.0e19
const K = 31.0613 / 4
const CUTOFF = 2.35
const NZ = 128
const N = 2
const J = 4

proton = BiKappa2(DENSITY, 5.5, 5.5, 1986.734, 993.367)
electron = Maxwellian(:e, DENSITY, 496.683)
data = gen_fv2d(proton)
in_plateau = abs.(sqrt(2) .* data.vz ./ data.vtz) .<= CUTOFF
fv = ifelse.(in_plateau, data.fv, zero(eltype(data.fv)))
alm = hermite_expansion(
    fv, data.vz, data.vx, data.vtz, data.vtx; Nz = NZ, Nx = 0,
).alm
proton_param = HHSolverParam(proton, B0; alm)
species = (proton_param, electron)
kx, kz = K .* sincos(THETA)
wci = abs(proton_param.wc)

function report_matrix(::Type{T}) where {T}
    return dispersion_matrix(
        species, T(B0), T(kx), T(kz), BOHH; N, J,
    )
end

values64 = eigvals(report_matrix(Float64))
values_double = eigvals(report_matrix(Double64))
values128 = setprecision(128) do
    eigvals(report_matrix(BigFloat))
end
reference = setprecision(192) do
    eigvals(report_matrix(BigFloat))
end

physical(values) = filter(values) do value
    abs(real(value) / wci) < 5 && abs(imag(value) / wci) < 2
end
reference_view = physical(reference)

function matched_roots(values)
    values_big = Complex{BigFloat}.(values)
    return map(reference_view) do target
        values_big[argmin(abs.(values_big .- target))]
    end
end

matched64 = matched_roots(values64)
matched_double = matched_roots(values_double)
matched128 = matched_roots(values128)
normalized(values) = values ./ wci

root_errors(values) = Float64.(log10.(abs.(values .- reference_view) ./ wci))
errors64 = root_errors(matched64)
errors_double = root_errors(matched_double)
errors128 = root_errors(matched128)
tick_indices = collect(1:8:length(reference_view))

function scatter_physical_roots!(axis)
    for (values, marker, color, label) in [
            (physical(reference), :circle, :black, "BigFloat(192) reference"),
            (physical(values64), :xcross, :crimson, "Float64"),
            (physical(values_double), :cross, :dodgerblue, "Double64"),
            (physical(values128), :star5, :seagreen, "BigFloat(128)"),
        ]
        roots = normalized(values)
        scatter!(
            axis, Float64.(real.(roots)), Float64.(imag.(roots));
            marker, markersize = 11, color, label
        )
    end
    return
end

figure = Figure(size = (1500, 470), fontsize = 17)
ax_roots = Axis(
    figure[1, 1],
    xlabel = "real(omega / Omega_p)",
    ylabel = "imag(omega / Omega_p)",
    title = "Physical-window solutions",
)
scatter_physical_roots!(ax_roots)
axislegend(ax_roots; position = :lb, labelsize = 12)

ax_zoom = Axis(
    figure[1, 2],
    xlabel = "real(omega / Omega_p)",
    ylabel = "imag(omega / Omega_p)",
    title = "Zoom: strongly shifted damped root",
)
scatter_physical_roots!(ax_zoom)
xlims!(ax_zoom, -0.3, 0.3)
ylims!(ax_zoom, -2.05, -0.3)
shift_index = argmax(abs.(matched64 .- reference_view))
shift_reference = normalized(reference_view[shift_index])
shift_float = normalized(matched64[shift_index])
lines!(
    ax_zoom, Float64.([real(shift_reference), real(shift_float)]),
    Float64.([imag(shift_reference), imag(shift_float)]);
    color = :crimson, linestyle = :dash, linewidth = 1.5
)
text!(ax_zoom, -0.29, -1.94; text = "reference", fontsize = 12, color = :black)
text!(ax_zoom, -0.12, -0.4; text = "Float64", fontsize = 12, color = :crimson)

ax_errors = Axis(
    figure[1, 3],
    xlabel = "reference root index",
    ylabel = "log10(abs(delta omega) / Omega_p)",
    title = "Matched-root difference",
    xticks = (tick_indices, string.(tick_indices)),
)
indices = collect(eachindex(reference_view))
scatterlines!(ax_errors, indices, errors64; marker = :xcross, color = :crimson, label = "Float64")
scatterlines!(ax_errors, indices, errors_double; marker = :cross, color = :dodgerblue, label = "Double64")
scatterlines!(ax_errors, indices, errors128; marker = :star5, color = :seagreen, label = "BigFloat(128)")
axislegend(ax_errors; position = :rb, labelsize = 12)

Label(
    figure[0, 1:3],
    "Sharp-edge firehose-like distribution: solution comparison (Nz = $NZ, N = $N, J = $J)",
    fontsize = 20,
    font = :bold,
)

output = pkgdir(PlasmaBO, "experiments", "arbitrary_precision", "arbitrary_precision_solution_comparison.svg")
save(output, figure)

vz = vec(data.vz[:, 1] ./ data.vtz)
vx = vec(data.vx[1, :] ./ data.vtx)
normalized_fv = fv ./ maximum(fv)
profile = vec(sum(fv; dims = 2))
profile ./= maximum(profile)
edge = CUTOFF / sqrt(2)

distribution = Figure(size = (1120, 430), fontsize = 17)
ax_distribution = Axis(
    distribution[1, 1],
    xlabel = "v_parallel / vtz",
    ylabel = "v_perp / vtp",
    title = "Hard-cutoff distribution",
)
heatmap!(
    ax_distribution, vz, vx,
    log10.(max.(normalized_fv, eps(Float64)));
    colormap = :viridis,
    colorrange = (-8, 0),
)
vlines!(ax_distribution, [-edge, edge]; color = :white, linestyle = :dash)

ax_profile = Axis(
    distribution[1, 2],
    xlabel = "v_parallel / vtz",
    ylabel = "integrated f / max",
    title = "Parallel lineout",
)
lines!(ax_profile, vz, profile; linewidth = 2, color = :dodgerblue)
vlines!(
    ax_profile, [-edge, edge]; color = :crimson, linestyle = :dash,
    label = "hard edge"
)
axislegend(ax_profile; position = :lb)
Label(
    distribution[0, 1:2],
    "Sharp-edge stress distribution used for solution comparison",
    fontsize = 20,
    font = :bold,
)
distribution_output = pkgdir(PlasmaBO, "experiments", "arbitrary_precision", "arbitrary_precision_sharp_edge_distribution.svg")
save(distribution_output, distribution)
err64 = Float64(abs(matched64[shift_index] - reference_view[shift_index]) / wci)
err_double = Float64(abs(matched_double[shift_index] - reference_view[shift_index]) / wci)
err128 = Float64(abs(matched128[shift_index] - reference_view[shift_index]) / wci)

results = joinpath(pkgdir(PlasmaBO), "experiments", "arbitrary_precision", "sharp_edge.toml")
open(results, "w") do io
    println(io, "# Auto-generated by precision_solution_plot.jl. Do not edit.")
    println(io, "reference = \"", string(shift_reference), "\"")
    println(io, "Float64_error = \"", Printf.@sprintf("%.4e", err64), "\"")
    println(io, "Double64_error = \"", Printf.@sprintf("%.4e", err_double), "\"")
    println(io, "BigFloat128_error = \"", Printf.@sprintf("%.4e", err128), "\"")
end
println("wrote $output")
println("wrote $distribution_output")
println("wrote $results")

end
