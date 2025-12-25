abstract type AbstractSolverParams end
abstract type AbstractSolverSpecies end

struct DispersionSolution{K, Θ, Ω}
    ks::K
    θs::Θ
    ωs::Ω
end

DispersionSolution(ks, ωs) = DispersionSolution(ks, nothing, ωs)

Base.eltype(::Type{DispersionSolution{K, Θ, Ω}}) where {K, Θ, Ω} = eltype(Ω)
Base.getindex(sol::DispersionSolution, i::Integer) = sol.ωs[i]
Base.getindex(sol::DispersionSolution, i::Integer, j::Integer) = sol.ωs[i, j]
Base.size(sol::DispersionSolution) = size(sol.ωs)
Base.axes(sol::DispersionSolution) = axes(sol.ωs)
Base.broadcastable(sol::DispersionSolution) = (sol,)
