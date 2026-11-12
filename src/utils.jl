# Valid until the enclosing `@no_escape` ends.
_temp_array(T, dims...) = isbitstype(T) ? Bumper.alloc!(Bumper.default_buffer(), T, dims...) : zeros(T, dims...)


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

Rotation sense of `v` about `B₀ ∥ +ẑ`: `:R` (electron-sense), `:L` (ion-sense) or `:linear`.
"""
function handedness(v, ω, k::Wavevector; threshold = 1.0e-5)
    E = @view v[(end - 5):(end - 3)]
    E⊥, Ey = _transverse_E(v, k), E[2]
    I = abs2(E⊥) + abs2(Ey)
    # Relative to |E|: v is unit-norm over all states, so its field part can be arbitrarily small.
    I <= threshold^2 * sum(abs2, E) && return :linear
    V = sign(real(ω)) * 2 * imag(conj(E⊥) * Ey) / I  # circular polarization degree
    return V > threshold ? :R : V < -threshold ? :L : :linear
end
