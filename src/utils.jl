macro temp_array(T, dims...)
    return esc(:(isbitstype($T) ? @alloc($T, $(dims...)) : zeros($T, $(dims...))))
end

_realtype(::Type{Complex{T}}) where {T} = T

# Execute function `f` with specified BLAS thread count, restoring previous setting afterward.
function solve_with_threads(f, nthreads)
    old = BLAS.get_num_threads()
    nthreads = min(old, nthreads)
    BLAS.set_num_threads(nthreads)
    return try
        f()
    finally
        BLAS.set_num_threads(old)
    end
end

function with_progress(f, prob; desc="Solving dispersion (k, θ)...", dt=1)
    θs = prob.θs
    ks = prob.ks
    carts = CartesianIndices((length(ks), length(θs)))
    return @showprogress dt = dt desc = desc Threads.@threads for id in carts
        ik, iθ = Tuple(id)
        k = ks[ik]
        θ = θs[iθ]
        kx, kz = k .* sincos(θ)
        f(ik, iθ, kx, kz)
    end
end

# Component of E perpendicular to k in the x-z plane, from an eigenvector ending in (E, B)
function _transverse_E(v, k::Wavevector)
    sinθ, cosθ = sincos(angle(k))
    return v[end - 5] * cosθ - v[end - 3] * sinθ
end

"""
    polarization_ratio(v, k::Wavevector)

Transverse polarization ratio `P = Ey / E⊥` of the eigenvector `v`, where `E⊥` is the transverse
electric field perpendicular to `k` in the x-z plane. `Ey / Ex` for parallel propagation.
"""
polarization_ratio(v, k::Wavevector) = v[end - 4] / _transverse_E(v, k)

"""
    handedness(v, ω, k::Wavevector; threshold = 1e-5)

Rotation sense in real time, about `B₀ ∥ +ẑ`, of the eigenvector `v` with complex frequency `ω`.

Returns `:R` (whistler-like, electron-sense), `:L` (Alfvén-like, ion-sense) or `:linear`.
"""
function handedness(v, ω, k::Wavevector; threshold = 1.0e-5)
    E_perp, Ey = _transverse_E(v, k), v[end - 4]
    abs(E_perp) < 1.0e-10 && abs(Ey) < 1.0e-10 && return :linear
    P = Ey / E_perp
    (abs(P) < threshold || 1 / abs(P) < threshold) && return :linear
    ellipticity = (real(ω) < 0 ? -1 : 1) * imag(P) / abs(P)
    return ellipticity > threshold ? :R : ellipticity < -threshold ? :L : :linear
end
