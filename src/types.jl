abstract type AbstractSolverParams end
abstract type AbstractSolverSpecies end
abstract type AbstractDispersionAlgorithm end


function _species(species)
    species isa AbstractVector && return species
    species isa Tuple && return species
    return (species,)
end

struct DispersionProblem{S, B, K}
    species::S
    B0::B
    kx::K
    kz::K

    function DispersionProblem(species, B0, kx, kz)
        species = _species(species)
        return new{typeof(species), typeof(B0), typeof(kx)}(species, B0, kx, kz)
    end
end

struct EnsembleProblem{S, B, K, Θ}
    species::S
    B0::B
    ks::K
    θs::Θ

    function EnsembleProblem(species, B0, ks, θs)
        species = _species(species)
        return new{typeof(species), typeof(B0), typeof(ks), typeof(θs)}(species, B0, ks, θs)
    end
end


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
