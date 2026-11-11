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

function _eigenroots(params, alg, ::Type{CT}) where {CT}
    n = matrix_size(alg, params)
    return function (k)
        ωs = @no_escape begin
            M = @temp_array(CT, n, n)
            dispersion_matrix!(zero!(M), params, k, alg)
            eigvals!(_for_eigvals(M))
        end
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
    params = prepare(alg, prob.plasma.populations, T)
    n = matrix_size(alg, params)
    M = dispersion_matrix!(zeros(Complex{T}, n, n), params, prob.k, alg)
    return M .*= u.omega
end
