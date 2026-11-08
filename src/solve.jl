include("solver/PBK.jl"); using .PBK: BOPBK
include("solver/HH.jl")
include("solver/fluid.jl")

zero!(M) = fill!(M, zero(eltype(M)))

# Promoted real (floating-point) element type carried by a problem.
function _realtype(::DispersionProblem{S, B, K}) where {S, B, K}
    return float(promote_type(B, K))
end
function _realtype(::EnsembleProblem{S, B, Ks, Θs}) where {S, B, Ks, Θs}
    return float(promote_type(B, eltype(Ks), eltype(Θs)))
end

_complex_type(pb) = Complex{_realtype(pb)}

_prepared(pb::DispersionProblem, alg) =
    DispersionProblem(prepare(alg, pb.species, pb.B0, _realtype(pb)), pb.B0, pb.kx, pb.kz)

function dispersion_matrix(pb::DispersionProblem, alg; kw...)
    pb = _prepared(pb, alg)
    n = matrix_size(alg, pb.species)
    return dispersion_matrix!(zeros(_complex_type(pb), n, n), pb, alg; kw...)
end

# LAPACK eigvals! accepts UnsafeArray views; GenericLinearAlgebra's Householder
# does not. Aliasing the Bumper buffer via `unsafe_wrap(Array, ...)` segfaults
# during a later GC mark (Householder intermediates outlive the @no_escape
# block), so we copy. The copy is negligible next to the extended-precision
# eigvals! itself.
_for_eigvals(M::AbstractMatrix{<:Union{ComplexF32, ComplexF64}}) = M
_for_eigvals(M::AbstractMatrix) = Matrix(M)

function solve(pb::DispersionProblem, alg)
    pb = _prepared(pb, alg)
    CT = _complex_type(pb)
    n = matrix_size(alg, pb.species)
    return @no_escape begin
        M = @temp_array(CT, n, n)
        dispersion_matrix!(zero!(M), pb, alg)
        eigvals!(_for_eigvals(M))
    end
end

function solve(pb::EnsembleProblem, alg)
    T = _realtype(pb)
    CT = Complex{T}
    ωs = Matrix{Vector{CT}}(undef, length(pb.ks), length(pb.θs))
    solve_with_threads(1) do
        species = prepare(alg, pb.species, pb.B0, T)
        n = matrix_size(alg, species)
        with_progress(pb) do ik, iθ, kx, kz
            @no_escape begin
                M = @temp_array(CT, n, n)
                sub = DispersionProblem(species, pb.B0, kx, kz)
                dispersion_matrix!(zero!(M), sub, alg)
                ωs[ik, iθ] = eigvals!(_for_eigvals(M))
            end
        end
    end
    return DispersionSolution(pb.ks, pb.θs, ωs)
end

"""
    solve(plasma::Plasma, k::Wavevector, alg) -> Vector of ω [rad/s]
    solve(plasma::Plasma, sweep::AngleSweep, alg) -> DispersionSolution

All modes at one wavevector [m⁻¹], or over `sweep.k` [m⁻¹] × `sweep.theta` [rad].

See also: [`BOHH`](@ref), [`BOPBK`](@ref), [`BOFluid`](@ref)
"""
solve(plasma::Plasma, k::Wavevector, alg::AbstractDispersionAlgorithm) = solve(DispersionProblem(plasma, k), alg)
solve(plasma::Plasma, sweep::AngleSweep, alg::AbstractDispersionAlgorithm) = solve(EnsembleProblem(plasma, sweep), alg)

"""
    dispersion_matrix(plasma::Plasma, k::Wavevector, alg)

The linear operator whose eigenvalues are the modes returned by [`solve`](@ref).
"""
dispersion_matrix(plasma::Plasma, k::Wavevector, alg::AbstractDispersionAlgorithm) =
    dispersion_matrix(DispersionProblem(plasma, k), alg)
