# Eigenvalue filtering and branch tracking for dispersion relations
using DataInterpolations

"""
    BranchPoint{T}
    BranchPoint(k, ω)

Initial point specification for dispersion branch tracking.

# Fields
- `k`: Wave vector at initial point
- `ω`: Frequency at initial point
- `locator`: Callable used to locate the seed eigenvalue at the initial point

# Example
```julia
# Track unstable branch starting at k=0.03 with γ=0.1721
point1 = BranchPoint(0.03, 0.1721im)     # Auto: track by Im

# Track real frequency mode at k=0.5 with ωᵣ=2.5
point2 = BranchPoint(0.5, 2.5)  # Auto: track by Re
```
"""
struct BranchPoint{T, W, L}
    k::T
    ω::W
    locator::L
end

abstract type AbstractSeedLocator end

struct SeedByReal <: AbstractSeedLocator end
struct SeedByImag <: AbstractSeedLocator end
struct SeedByAbs <: AbstractSeedLocator end

(::SeedByReal)(xs, x0) = argmin(x -> abs(real(x) - real(x0)), xs)
(::SeedByImag)(xs, x0) = argmin(x -> abs(imag(x) - imag(x0)), xs)
(::SeedByAbs)(xs, x0) = argmin(x -> abs(x - x0), xs)

_default_seed_locator(ω) =
    iszero(imag(ω)) ? SeedByReal() : (iszero(real(ω)) ? SeedByImag() : SeedByAbs())

BranchPoint(k, ω) = BranchPoint(k, ω, _default_seed_locator(ω))

struct SurfaceBranchPoint{T, W, L}
    k::T
    θ::T
    ω::W
    locator::L
end

SurfaceBranchPoint(k, θ, ω) =
    SurfaceBranchPoint(k, θ, ω, _default_seed_locator(ω))

# Interpolate complex values by treating real and imaginary parts independently
function interpolate_complex(z_prev, x_prev, x_new; extrapolation = ExtrapolationType.Extension)
    return if length(z_prev) == 1
        z_prev[1]
    elseif length(z_prev) == 2
        # Linear interpolation
        z1, z2 = z_prev
        x1, x2 = x_prev
        z1 + (z2 - z1) * (x_new - x1) / (x2 - x1)
    else
        # PCHIP interpolation
        itp_real = PCHIPInterpolation(real.(z_prev), x_prev; extrapolation)
        itp_imag = PCHIPInterpolation(imag.(z_prev), x_prev; extrapolation)
        complex(itp_real(x_new), itp_imag(x_new))
    end
end

track(sol, point::NTuple{2, Any}) =
    track(sol, BranchPoint(point...))

track(sol, point::NTuple{3, Any}) =
    track(sol, SurfaceBranchPoint(point...))

"""
    track(solution, point)

Track a single dispersion branch across parameter space from initial `point`.

The algorithm:
1. Start at the k-point nearest to the initial point
2. Track bidirectionally using interpolation to predict ω at each k
3. Select nearest eigenvalue to prediction

# Example
```julia
initial = BranchPoint(0.03, 0.1721im)
k_branch, ω_branch = track(solution, initial)
```

See also: [`BranchPoint`](@ref)
"""
track(sol, point::BranchPoint) = _track(sol.ks, sol, point)


function _track(ks, result, k, ω_seed)
    nk = length(ks)
    start_idx = argmin(abs.(ks .- k))
    ω_branch = similar(eltype(result), nk)
    ω_branch[start_idx] = ω_seed

    # Track in both directions
    track_from_index!(ω_branch, ks, result, start_idx, (start_idx + 1):nk)
    track_from_index!(ω_branch, ks, result, start_idx, (start_idx - 1):-1:1)

    return ks, ω_branch
end

function _track(ks, result, point::BranchPoint, ω_seed = nothing)
    # Initialize at starting point
    start_idx = argmin(abs.(ks .- point.k))
    ω_seed = @something ω_seed point.locator(result[start_idx], point.ω)
    return _track(ks, result, ks[start_idx], ω_seed)
end


# Track along one direction from a starting point
function track_from_index!(ω_branch, ks, results, start_idx, indices)
    isempty(indices) && return

    tracked = Int[start_idx]
    is_forward = first(indices) > start_idx
    for i in indices
        # Get previously tracked indices in this direction
        prev = is_forward ? filter(<(i), tracked) : filter(>(i), tracked)
        isempty(prev) && continue

        # Predict ω at new k using interpolation
        sort!(prev)
        ω_pred = interpolate_complex(ω_branch[prev], ks[prev], ks[i])

        # Find nearest eigenvalue to prediction
        ω_branch[i] = argmin(x -> abs(x - ω_pred), results[i])
        push!(tracked, i)
    end
    return
end

function track(sol, point::SurfaceBranchPoint)
    ks, θs = sol.ks, sol.θs
    nk, nθ = length(ks), length(θs)

    # Find starting indices
    idx_k0 = argmin(abs.(ks .- point.k))
    idx_θ0 = argmin(abs.(θs .- point.θ))

    ω_surface = Matrix{ComplexF64}(undef, nk, nθ)

    ω_seed = point.locator(sol.ωs[idx_k0, idx_θ0], point.ω)
    # Fill initial θ slice
    fill_θ_slice!(ω_surface, sol, ks, idx_θ0, point.k, ω_seed)

    # Track in both θ directions
    track_θ_direction!(ω_surface, sol, ks, θs, idx_k0, idx_θ0, (idx_θ0 + 1):nθ, point)
    track_θ_direction!(ω_surface, sol, ks, θs, idx_k0, idx_θ0, (idx_θ0 - 1):-1:1, point)

    return ks, θs, ω_surface
end

function fill_θ_slice!(ω_surface, result, ks, iθ, k, ω_seed)
    ωsθ = result.ωs[:, iθ]
    _, ωb = _track(ks, ωsθ, k, ω_seed)
    return ω_surface[:, iθ] .= ωb
end

function track_θ_direction!(ω_surface, result, ks, θs, idx_k0, idx_θ0, indices, point)
    isempty(indices) && return

    tracked = Int[idx_θ0]
    is_forward = first(indices) > idx_θ0

    for iθ in indices
        # Get previously tracked indices
        prev = is_forward ? filter(<(iθ), tracked) : filter(>(iθ), tracked)
        isempty(prev) && continue

        # Predict ω at k0 for this θ
        sort!(prev)
        ω_at_k0 = view(ω_surface, idx_k0, :)
        ω_pred = interpolate_complex(ω_at_k0[prev], θs[prev], θs[iθ])
        ω_seed = argmin(x -> abs(x - ω_pred), result.ωs[idx_k0, iθ])
        fill_θ_slice!(ω_surface, result, ks, iθ, point.k, ω_seed)
        push!(tracked, iθ)
    end
    return
end
