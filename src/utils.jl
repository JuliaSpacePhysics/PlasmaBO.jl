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
relative to the background magnetic field direction (which is along the +z direction).

Returns:
* `:R` (Right-handed/whistler-like/electron-sense)
* `:L` (Left-handed/Alfvén-like/ion-sense)
* `:linear` (Linearly polarized within the threshold)
"""
function handedness(v; threshold=1e-5)
    P = polarization_ratio(v)
    im_P = imag(P)
    if im_P > threshold
        return :R
    elseif im_P < -threshold
        return :L
    else
        return :linear
    end
end

