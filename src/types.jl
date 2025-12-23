abstract type AbstractSolverParams end
abstract type AbstractSolverSpecies end

struct DispersionSolution{K, Ω}
    ks::K
    ωs::Ω
end
