abstract type AbstractDispersionAlgorithm end

"""
    DispersionProblem(plasma::Plasma, k::Wavevector)

All modes at one wavevector [m⁻¹], `B₀ ∥ ẑ`.
"""
struct DispersionProblem{S, B, K}
    species::S
    B0::B
    kx::K
    kz::K

    # `species`: `Population`s, or an algorithm's `prepare`d parameters
    function DispersionProblem(species, B0, kx::K1, kz::K2) where {K1, K2}
        K = promote_type(K1, K2)
        return new{typeof(species), typeof(B0), K}(species, B0, K(kx), K(kz))
    end
end

DispersionProblem(plasma::Plasma, k::Wavevector) = DispersionProblem(plasma.populations, norm(plasma.B0), k.k_perp, k.k_para)

# All modes over wavenumbers [m⁻¹] × propagation angles [rad] from `B₀`.
struct EnsembleProblem{S, B, K, Θ}
    species::S
    B0::B
    ks::K
    θs::Θ
end

EnsembleProblem(plasma::Plasma, g::AngleSweep) =
    EnsembleProblem(plasma.populations, norm(plasma.B0), g.k, g.theta)

struct DispersionSolution{K, Θ, Ω}
    k::K
    theta::Θ
    omega::Ω
end

DispersionSolution(k, omega) = DispersionSolution(k, nothing, omega)

Base.eltype(::Type{DispersionSolution{K, Θ, Ω}}) where {K, Θ, Ω} = eltype(Ω)
Base.getindex(sol::DispersionSolution, i::Integer) = sol.omega[i]
Base.getindex(sol::DispersionSolution, i::Integer, j::Integer) = sol.omega[i, j]
Base.size(sol::DispersionSolution) = size(sol.omega)
Base.axes(sol::DispersionSolution) = axes(sol.omega)
Base.broadcastable(sol::DispersionSolution) = (sol,)
