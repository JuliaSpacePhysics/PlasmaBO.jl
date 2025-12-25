module PlasmaBOMakieExt

using PlasmaBO
using PlasmaBO: DispersionSolution
using Makie
import PlasmaBO: plot_branches

struct FigureAxes{A}
    figure::Figure
    axes::A
end

Base.iterate(fg::FigureAxes, args...) = iterate((fg.figure, fg.axes), args...)
Base.show(io::IO, m::MIME, fg::FigureAxes) = show(io, m, fg.figure)
Base.show(io::IO, ::MIME"text/plain", fg::FigureAxes) = print(io, "FigureAxes()")
Base.showable(m::MIME, fg::FigureAxes) = showable(m, fg.figure)

"""
    plot_dispersion(solution; kw...)

Plot a dispersion solution showing eigenvalues across parameter space.

For 1D solutions (k only), creates scatter plots of Re(ω) and Im(ω) vs k.
For 2D solutions (k, θ), creates surface plots.

# Keyword Arguments
- `normalize`: Optional named tuple `(k=kn, ω=ωn)` for normalization
- Additional keyword arguments are passed to the plotting functions
"""
function Makie.plot(sol::DispersionSolution, args...; kwargs...)
    pf = isnothing(sol.θs) ? plot_1d : plot_2d
    return pf(sol, args...; kwargs...)
end

const XLABEL = L"k/k_n"
const LABEL_REAL = L"ω_r \; / \; ω_n"
const LABEL_IMAG = L"ω_i \; / \; ω_n"

function plot_1d(solution, kn, ωn; color = Cycled(1), xlabel = XLABEL, kwargs...)
    return with_2d_axes(; xlabel) do ax1, ax2
        for (k, ωs) in zip(solution.ks, solution.ωs)
            k_norm = fill(k / kn, length(ωs))
            scatter!(ax1, k_norm, real.(ωs) ./ ωn; markersize = 5, color = color, kwargs...)
            scatter!(ax2, k_norm, imag.(ωs) ./ ωn; markersize = 5, color = color, kwargs...)
        end
    end
end

function with_2d_axes(func!, args...; xlabel = XLABEL, figure = (;), by_col = false)
    layout = by_col ? ((1, 1), (1, 2)) : ((1, 1), (2, 1))
    fig = Figure(; figure...)
    ax1 = Axis(fig[layout[1]...]; xlabel, ylabel = LABEL_REAL)
    ax2 = Axis(fig[layout[2]...]; xlabel, ylabel = LABEL_IMAG)
    func!(ax1, ax2, args...)
    !by_col && hidexdecorations!(ax1, grid = false)
    return FigureAxes(fig, (ax1, ax2))
end

function with_3d_axes(func!, args...; xlabel = XLABEL, ylabel = "θ (rad)", figure = (;), by_col = false)
    layout = by_col ? ((1, 1), (1, 2)) : ((1, 1), (2, 1))
    fig = Figure(; figure...)
    ax1 = Axis3(fig[layout[1]...]; xlabel, ylabel, zlabel = LABEL_REAL)
    ax2 = Axis3(fig[layout[2]...]; xlabel, ylabel, zlabel = LABEL_IMAG)
    func!(ax1, ax2, args...)
    !by_col && hidexdecorations!(ax1, grid = false)
    return FigureAxes(fig, (ax1, ax2))
end

function plot_2d(sol, kn, ωn; xlabel = XLABEL, color = Cycled(1), figure = (;), kwargs...)
    return with_3d_axes(; figure, xlabel) do ax1, ax2
        for i in axes(sol.ωs, 1), j in axes(sol.ωs, 2)
            ωs = sol.ωs[i, j]
            k_norm = fill(sol.ks[i] / kn, length(ωs))
            θ_scaled = fill(sol.θs[j], length(ωs))

            # Extract real and imaginary parts for this point
            z_real = real(ωs) / ωn
            z_imag = imag(ωs) / ωn

            scatter!(ax1, k_norm, θ_scaled, z_real; markersize = 5, color, kwargs...)
            scatter!(ax2, k_norm, θ_scaled, z_imag; markersize = 5, color, kwargs...)
        end
    end
end

"""
    plot_branches(branches, kn, ωn; kw...)

Plot dispersion `branches` from branch tracking.
"""

function PlasmaBO.plot_branches(branches, args...; kwargs...)
    f = length(first(branches)) == 2 ? _plot_branches_1d : _plot_branches_2d
    return f(branches, args...; kwargs...)
end

function _plot_branches_1d(branches, kn, ωn; figure = (;), xlabel = XLABEL, kwargs...)
    return with_2d_axes(; figure, xlabel) do ax1, ax2
        for (i, (k, ω)) in enumerate(branches)
            lines!(ax1, k ./ kn, real.(ω) ./ ωn; label = "Branch $i", kwargs...)
            lines!(ax2, k ./ kn, imag.(ω) ./ ωn; kwargs...)
        end
    end
end

function _plot_branches_2d(branches, kn, ωn; figure = (;), kwargs...)
    return with_3d_axes(; figure) do ax1, ax2
        for (_, (k, θ, ω)) in enumerate(branches)
            surface!(ax1, k ./ kn, θ, real.(ω) ./ ωn; kwargs...)
            # wireframe!(ax1, k ./ kn, θ, real.(ω) ./ ωn; color = :black, kwargs...)
            surface!(ax2, k ./ kn, θ, imag.(ω) ./ ωn; kwargs...)
            # wireframe!(ax2, k ./ kn, θ, imag.(ω) ./ ωn; color = :black, kwargs...)
        end
    end
end

end # module
