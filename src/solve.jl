abstract type AbstractDispersionAlgorithm end

zero!(M) = fill!(M, zero(eltype(M)))

# The 3×3 coupling of a velocity triple (x, y, z) to E.
function _coupling!(C, (p11, p12, p13, p22, p23, p33))
    C[1, 1], C[1, 2], C[1, 3] = p11, p12, p13
    C[2, 1], C[2, 2], C[2, 3] = -p12, p22, p23
    C[3, 1], C[3, 2], C[3, 3] = p13, -p23, p33
    return C
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

# GenericLinearAlgebra's eigvals! fails on Bumper's UnsafeArray, and an `unsafe_wrap(Array, ...)`
# alias crashes GC (its intermediates outlive @no_escape), so non-LAPACK types get a copy.
_for_eigvals(M::AbstractMatrix{<:Union{ComplexF32, ComplexF64}}) = M
_for_eigvals(M::AbstractMatrix) = Matrix(M)

function pointsolver(prob, alg::AbstractDispersionAlgorithm)
    T = realtype(prob)
    params = prepare(alg, prob.plasma.populations, T)
    return function (k)
        ωs = @no_escape eigvals!(_for_eigvals(_operator(_temp_array, params, k, alg, Complex{T})))
        return (; omega = ωs, resid = nothing, nevals = 1, converged = true)
    end
end

function _operator(alloc, params, k, alg, ::Type{CT}) where {CT}
    n = matrix_size(alg, params)
    return dispersion_matrix!(zero!(alloc(CT, n, n)), params, k, alg)
end

"""
    dispersion_matrix(plasma, k::Wavevector, alg; ref)

The linear operator whose eigenvalues are the roots `ω`, in the units of `plasma` (see
`DispersionProblem` for `ref`). Its eigenvectors end in the fields `(E, B)` with `c = 1`.
"""
function dispersion_matrix(plasma, k::Wavevector, alg::AbstractDispersionAlgorithm; ref...)
    prob, u = normalized(DispersionProblem(plasma, k; ref...))
    T = realtype(prob)
    M = _operator(zeros, prepare(alg, prob.plasma.populations, T), prob.k, alg, Complex{T})
    return M .*= u.omega
end
