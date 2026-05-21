# Complementary error function via libopenlibm — same path SpecialFunctions uses for Float64.
erfc(x::Float64) = ccall((:erfc, Base.Math.libm), Float64, (Float64,), x)
erfc(x::Real) = erfc(float(Float64(x)))

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

function with_progress(f, prob; desc = "Solving dispersion (k, θ)...", dt = 1)
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

"""
    electric_field(eigenvector)

Extract the electric field components `(Ex, Ey, Ez)` from an eigenvector.
"""
function electric_field(v)
    n = length(v)
    return (v[n-5], v[n-4], v[n-3])
end

"""
    magnetic_field(eigenvector)

Extract the magnetic field components `(Bx, By, Bz)` from an eigenvector.
"""
function magnetic_field(v)
    n = length(v)
    return (v[n-2], v[n-1], v[n])
end

"""
    polarization_ratio(v)

Return the transverse polarization ratio `P = Ey / Ex` of the eigenvector.
"""
function polarization_ratio(v)
    Ex, Ey, _ = electric_field(v)
    return Ey / Ex
end

"""
    handedness(v; threshold=1e-5)

Determine the polarization handedness of the wave mode represented by the eigenvector `v`
relative to the background magnetic field direction (along the +z direction), assuming a 
positive frequency (ω > 0) wave.
"""
function handedness(v; threshold=1e-5)
    return handedness(v, 1.0; threshold=threshold)
end

"""
    handedness(v, ω; threshold=1e-5)

Determine the physical polarization handedness of the wave mode represented by the eigenvector `v`
with complex frequency `ω` relative to the background magnetic field direction (+z direction).

We assume the wave phase convention is e^{i(k_z z - ω t)} with wave propagation along +z.
This physical handedness correctly factors in the sign of real(ω) to determine the real-time 
rotation sense of the physical field.

Returns:
* `:R` (Right-handed/whistler-like/electron-sense)
* `:L` (Left-handed/Alfvén-like/ion-sense)
* `:linear` (Linearly polarized or near-linear within the threshold)
"""
function handedness(v, ω; threshold=1e-5)
    Ex, Ey, _ = electric_field(v)
    abs_Ex = abs(Ex)
    abs_Ey = abs(Ey)
    if abs_Ex < 1e-10 && abs_Ey < 1e-10
        return :linear
    end
    P = Ey / Ex
    abs_P = abs(P)
    if abs_P < threshold || 1/abs_P < threshold
        return :linear
    end
    # Physical handedness factors in the sign of real(ω) to reflect true time-domain rotation
    sgn = sign(real(ω))
    # If ω is exactly 0 (or pure imaginary), default to positive frequency convention
    sgn = sgn == 0 ? 1.0 : sgn
    ellipticity = (sgn * imag(P)) / abs_P
    if ellipticity > threshold
        return :R
    elseif ellipticity < -threshold
        return :L
    else
        return :linear
    end
end

