abstract type AbstractDispersionAlgorithm end

zero!(M) = fill!(M, zero(eltype(M)))

# Rows r, r + SNJ1, r + 2SNJ1 couple to E through [p11 p12 p13; -p12 p22 p23; p13 -p23 p33].
function _couple_E!(M, r, SNJ1, SNJ3, (p11, p12, p13, p22, p23, p33))
    M[r, (SNJ3 + 1):(SNJ3 + 3)] .= (p11, p12, p13)
    M[r + SNJ1, (SNJ3 + 1):(SNJ3 + 3)] .= (-p12, p22, p23)
    M[r + 2SNJ1, (SNJ3 + 1):(SNJ3 + 3)] .= (p13, -p23, p33)
    return M
end

# Maxwell's curl equations in units with c = 1
function _B_E_part!(M, idx, kx, kz)
    # E(B) coupling: Ampère's law
    M[idx + 1, idx + 5] = kz
    M[idx + 2, idx + 4] = -kz
    M[idx + 2, idx + 6] = kx
    M[idx + 3, idx + 5] = -kx

    # B(E) coupling: Faraday's law
    M[idx + 4, idx + 2] = -kz
    M[idx + 5, idx + 1] = kz
    M[idx + 5, idx + 3] = -kx
    M[idx + 6, idx + 2] = kx
    return M
end

include("solver/PBK.jl")
include("solver/HH.jl")
include("solver/fluid.jl")

# LAPACK eigvals! accepts UnsafeArray views; GenericLinearAlgebra's Householder
# does not. Aliasing the Bumper buffer via `unsafe_wrap(Array, ...)` segfaults
# during a later GC mark (Householder intermediates outlive the @no_escape
# block), so we copy. The copy is negligible next to the extended-precision
# eigvals! itself.
_for_eigvals(M::AbstractMatrix{<:Union{ComplexF32, ComplexF64}}) = M
_for_eigvals(M::AbstractMatrix) = Matrix(M)

# The static parameter keeps the closure's element type inferable.
pointsolver(prob, alg::AbstractDispersionAlgorithm) =
    _eigenroots(prepare(alg, prob.plasma.populations, realtype(prob)), alg, Complex{realtype(prob)})

_bumper_alloc(T, dims...) = isbitstype(T) ? Bumper.alloc!(Bumper.default_buffer(), T, dims...) : Array{T}(undef, dims...)

function _operator(alloc, params, k, alg, ::Type{CT}) where {CT}
    n = matrix_size(alg, params)
    return dispersion_matrix!(zero!(alloc(CT, n, n)), params, k, alg)
end

function _eigenroots(params, alg, ::Type{CT}) where {CT}
    return function (k)
        ωs = @no_escape eigvals!(_for_eigvals(_operator(_bumper_alloc, params, k, alg, CT)))
        return (; omega = ωs, resid = nothing, nevals = 1, converged = true)
    end
end

"""
    dispersion_matrix(plasma, k::Wavevector, alg; ref)

The linear operator whose eigenvalues are the roots `ω`, in the units of `plasma` (see
`DispersionProblem` for `ref`). Its eigenvectors end in the fields `(E, B)` with `c = 1`.
"""
function dispersion_matrix(plasma, k::Wavevector, alg::AbstractDispersionAlgorithm; ref...)
    prob, u = normalized(DispersionProblem(plasma, k; ref...))
    T = realtype(prob)
    M = _operator((S, dims...) -> zeros(S, dims...), prepare(alg, prob.plasma.populations, T), prob.k, alg, Complex{T})
    return M .*= u.omega
end

# Triples r = 1:m: rows r, r + m, r + 2m hold (x, y, z) of ẋᵣ = ξᵣ xᵣ + Bᵣ E; then E (with Ė = -Σᵣ xᵣ + k × B) and B.
function _full_operator(alloc, ξ, B, kx, kz)
    m = length(ξ)
    M = zero!(alloc(eltype(ξ), 3m + 6, 3m + 6))
    for r in 1:m, c in 0:2
        M[r + c * m, r + c * m] = ξ[r]
        M[r + c * m, (3m + 1):(3m + 3)] .= view(B, c + 1, :, r)
        M[3m + 1 + c, r + c * m] = -1
    end
    return _B_E_part!(M, 3m, kx, kz)
end
