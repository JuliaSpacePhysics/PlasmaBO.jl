abstract type AbstractSolverParams end
abstract type AbstractSolverSpecies end

struct DispersionSolution{K, Ω}
    ks::K
    ωs::Ω
end

# Distribution function generation and Hermite expansion utilities
"""
    Species

Species parameters for plasma dispersion relation solver.

# Fields
- `q`: Charge in units of electron charge (e.g., -1 for electron, 1 for proton)
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity (in units of c)
- `vdr`: Perpendicular ring beam drift velocity (in units of c)
- `aslm`: Hermite expansion coefficients (optional, computed if not provided)
"""
struct Species{T}
    q::T           # Charge (units of e)
    m::T           # Mass (units of mp)
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    vdr::T         # Perpendicular drift velocity (c)
    aslm::Matrix{T}  # Hermite expansion coefficients
end

"""
    Species(q, m, n, Tz, Tp; vdz=0.0, vdr=0.0, distribution=:maxwellian, kwargs...)

Create a Species with automatic Hermite coefficient computation.

# Distributions
- `:maxwellian` or `:bi_maxwellian`: Drift bi-Maxwellian ring beam
- `:bi_kappa`: Bi-kappa distribution
- `:product_bi_kappa`: Product bi-kappa distribution
- `:kappa_maxwellian`: Kappa-Maxwellian distribution
- `:shell`: Shell distribution
- `:slowing_down`: Slowing down distribution

# Keyword Arguments
- `kappa`: Kappa index for kappa distributions (default: 5.5)
- `Nz`, `Nx`: Maximum Hermite indices (default: 16)

See also: [`expand_fv2d`](@ref), [`gen_fv2d`](@ref)
"""
function Species(q, m, n, Tz, Tp = Tz; vdz = 0.0, vdr = 0.0, distribution = :maxwellian, kwargs...)
    aslm = if distribution == :maxwellian || distribution == :bi_maxwellian
        # For Maxwellian, aslm is simply 1.0
        ones(eltype(q), 1, 1)
    else
        # Generate distribution on grid and expand in Hermite basis
        hermite_coefficients(Tz, Tp, vdz, vdr, m; distribution, kwargs...)
    end

    return Species(promote(q, m, n, Tz, Tp, vdz, vdr)..., aslm)
end
