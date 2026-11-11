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

# The operator with auxiliary triples r (ẋᵣ = ξᵣ xᵣ + Bᵣ E, Ė = -Σᵣ xᵣ + k × B) reduced to the states
# the fields couple to: at each distinct pole, as many as the rank of the summed residue, with pivots
# below the poles' `tol` cut (Tian, Xie & Wei 2026, arXiv:2610.41450, §3.1). Roots are unchanged
# (their eq. 15); the dropped eigenvalues sit exactly on the poles.
function _reduced_operator(alloc, ξ, B, tol, kx, kz)
    CT = eltype(ξ)
    m = length(ξ)
    # + 0 maps -0.0 to 0.0, so equal poles sort into contiguous runs.
    perm = sortperm(ξ; by = z -> (real(z) + 0, imag(z) + 0))
    ξr = alloc(CT, 3m)
    L = alloc(CT, 3, 3m)    # E rows ← reduced states
    Rf = alloc(CT, 3m, 3)   # reduced states ← E
    Q = alloc(CT, 3, 3)
    r = 0
    i = 1
    while i <= m
        pole = ξ[perm[i]]
        fill!(Q, zero(CT))
        cut = zero(real(CT))
        while i <= m && ξ[perm[i]] == pole
            Q .-= view(B, :, :, perm[i])
            cut = max(cut, tol[perm[i]])
            i += 1
        end
        # Complete-pivoting elimination of Q: each pivot Q[a, b] peels the rank-1 term Q[:, b] Q[a, :] / Q[a, b].
        for _ in 1:3
            p, idx = findmax(abs, Q)
            p > cut || break
            a, b = Tuple(idx)
            s = sqrt(Q[a, b])   # split the pivot's size evenly between the two couplings
            r += 1
            ξr[r] = pole
            for c in 1:3
                L[c, r] = Q[c, b] / s
                Rf[r, c] = Q[a, c] / s
            end
            for c in 1:3, d in 1:3
                Q[d, c] -= L[d, r] * Rf[r, c]
            end
        end
    end
    M = zero!(alloc(CT, r + 6, r + 6))
    for t in 1:r
        M[t, t] = ξr[t]
    end
    M[1:r, (r + 1):(r + 3)] .= view(Rf, 1:r, :)
    M[(r + 1):(r + 3), 1:r] .= view(L, :, 1:r)
    return _B_E_part!(M, r, kx, kz)
end
