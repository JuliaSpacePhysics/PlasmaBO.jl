# Eigenvalue filtering and branch tracking for dispersion relations
using DataInterpolations

"""
    BranchPoint{T}
    BranchPoint(k, ω)

Initial point specification for dispersion branch tracking.

# Fields
- `k`: Wave vector at initial point
- `ω`: Frequency at initial point

Optional field:
- `track_by_real`: If true, use Re(ω) for initial point search; if false, use Im(ω), default is true if Re(ω) != 0.

# Example
```julia
# Track unstable branch starting at k=0.03 with γ=0.1721
point1 = BranchPoint(0.03, 0.1721im)     # Auto: track by Im

# Track real frequency mode at k=0.5 with ωᵣ=2.5
point2 = BranchPoint(0.5, 2.5)  # Auto: track by Re
```
"""
struct BranchPoint{T}
    k::T
    ω::Complex{T}
    track_by_real::Bool
end

_k(b::BranchPoint) = b.k
_ω(b::BranchPoint) = b.ω
_track_by_real(b::BranchPoint) = b.track_by_real
_k(tuple) = tuple[1]
_ω(tuple) = tuple[2]
_track_by_real(tuple) = real(tuple[2]) != 0


function BranchPoint(k::T, ω::Complex{T}) where {T}
    return BranchPoint{T}(k, ω, real(ω) != 0)
end

BranchPoint(k, ω::Real) = BranchPoint(k, complex(ω))

"""
    track_dispersion_branch(ks, ωs, point)

Track a single dispersion branch across parameter space from initial `point` using interpolation given wave vectors `ks` and a collection of eigenvalues `ωs`.

Returns parameter values and eigenvalues along the branch

The branch tracking algorithm is implemented as follows:
1. Start with initial point (k₀, ω₀)
2. Use interpolation to predict ω at next k
3. Find nearest eigenvalue to prediction
4. Continue building the branch

# Example
```julia
initial = BranchPoint(0.03, 0.1721im)
k_branch, ω_branch = track_dispersion_branch(k_scan, ω_all, initial)
```

See also: [`BranchPoint`](@ref), [`track_dispersion_branches`](@ref)
"""
function track_dispersion_branch(ks, ωs, point)
    nk = length(ks)
    T = eltype(ks)
    k = _k(point)
    ω = _ω(point)
    # Initialize output arrays
    k_branch = zeros(T, nk)
    ω_branch = zeros(Complex{T}, nk)

    # Get eigenvalues as matrix or compute on-the-fly
    get_eigenvalues = if ωs isa Function
        # Compute eigenvalues on demand
        ωs
    else
        # Pre-computed eigenvalues
        i -> ωs[i]
    end

    # Find starting index closest to initial k
    start_idx = argmin(abs.(ks .- k))
    k_start = ks[start_idx]

    # Find starting eigenvalue
    ω_at_start = get_eigenvalues(start_idx)
    idx = if _track_by_real(point)
        argmin(abs.(real.(ω_at_start) .- real(ω)))
    else
        argmin(abs.(imag.(ω_at_start) .- imag(ω)))
    end
    ω_start = ω_at_start[idx]

    # Store starting point
    k_branch[start_idx] = k_start
    ω_branch[start_idx] = ω_start

    # Track in both directions from starting point
    tracked_indices = Set([start_idx])

    # Helper function to track in one direction
    function track_direction(indices)
        # Track indices for this direction only (like MATLAB's jjpa)
        direction_indices = Int[start_idx]

        for i in indices
            if i in tracked_indices
                continue
            end

            # Get previous tracked points from THIS direction only
            # Forward direction: indices[1] < indices[end], want j < i
            # Backward direction: indices[1] > indices[end], want j > i
            is_forward = indices[1] < indices[end]
            prev_indices = sort(
                [
                    j for j in direction_indices if
                        (is_forward ? j < i : j > i)
                ]
            )

            if length(prev_indices) == 0
                continue
            end

            # Predict ω at current k using interpolation
            if length(prev_indices) == 1
                # Nearest neighbor for first step
                ω_pred = ω_branch[prev_indices[1]]
            elseif length(prev_indices) == 2
                # Linear extrapolation
                k1, k2 = k_branch[prev_indices[1]], k_branch[prev_indices[2]]
                ω1, ω2 = ω_branch[prev_indices[1]], ω_branch[prev_indices[2]]
                k_new = ks[i]
                ω_pred = ω1 + (ω2 - ω1) * (k_new - k1) / (k2 - k1)
            else
                # PCHIP interpolation requires at least 3 points
                k_prev = k_branch[prev_indices]
                ω_prev = ω_branch[prev_indices]
                ω_pred = z_interpolate(ω_prev, k_prev, ks[i])
            end

            # Find nearest eigenvalue to prediction
            ω_at_i = get_eigenvalues(i)
            idx = argmin(abs.(ω_at_i .- ω_pred))

            # Store result
            k_branch[i] = ks[i]
            ω_branch[i] = ω_at_i[idx]
            push!(tracked_indices, i)
            push!(direction_indices, i)  # Add to direction-specific list
        end
        return
    end

    # Track to the right and left
    track_direction((start_idx + 1):nk)
    track_direction((start_idx - 1):-1:1)

    return k_branch, ω_branch
end

# It seems that we need to interpolate real and imaginary parts independently
function z_interpolate(z_prev, t, t_new; extrapolation = ExtrapolationType.Extension)
    T = PCHIPInterpolation
    itp_real = T(real.(z_prev), t; extrapolation)
    itp_imag = T(imag.(z_prev), t; extrapolation)
    return complex(itp_real(t_new), itp_imag(t_new))
end

"""
    track_dispersion_branches(ks, eigenvalues, points)

Track multiple dispersion branches from initial `points` simultaneously given `ks` and `eigenvalues`.

Returns a vector of (k_branch, ω_branch) tuples, one per branch

See also: [`track_dispersion_branch`](@ref), [`BranchPoint`](@ref)
"""
function track_dispersion_branches(ks, eigenvalues, points)
    return map(points) do point
        track_dispersion_branch(ks, eigenvalues, point)
    end
end


track_dispersion_branches(result::DispersionSolution, points) =
    track_dispersion_branches(result.ks, result.ωs, points)
