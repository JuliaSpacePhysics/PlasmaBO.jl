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

# Emits ProgressLogging-protocol messages when `progress=true`.
function with_progress(f, prob;
        progress = false,
        name = "Solving dispersion (k, θ)...",
        progress_steps = 1)
    θs = prob.θs
    ks = prob.ks
    carts = CartesianIndices((length(ks), length(θs)))
    if !progress
        Threads.@threads for id in carts
            ik, iθ = Tuple(id)
            kx, kz = ks[ik] .* sincos(θs[iθ])
            f(ik, iθ, kx, kz)
        end
        return
    end
    total = length(carts)
    counter = Threads.Atomic{Int}(0)
    @withprogress name = name begin
        Threads.@threads for id in carts
            ik, iθ = Tuple(id)
            kx, kz = ks[ik] .* sincos(θs[iθ])
            f(ik, iθ, kx, kz)
            done = Threads.atomic_add!(counter, 1) + 1
            if done == total || done % progress_steps == 0
                @logprogress done / total
            end
        end
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
    polarization_ratio(v, kx, kz)
    polarization_ratio(v; kx=0.0, kz=0.0)

Return the transverse polarization ratio `P = Ey / E_perp` of the eigenvector `v`, where `E_perp`
is the transverse electric field component perpendicular to the wavevector `k` in the `x-z` plane:
`E_perp = Ex * cos(θ) - Ez * sin(θ)`.

If both `kx` and `kz` are zero (or omitted), it defaults to parallel propagation (`kx = 0.0, kz = 1.0`),
where the ratio simplifies to `Ey / Ex`.
"""
function polarization_ratio(v, kx, kz)
    Ex, Ey, Ez = electric_field(v)
    k = sqrt(kx^2 + kz^2)
    if k == 0.0
        return Ey / Ex
    else
        cos_θ = kz / k
        sin_θ = kx / k
        E_perp = Ex * cos_θ - Ez * sin_θ
        return Ey / E_perp
    end
end

function polarization_ratio(v; kx=0.0, kz=0.0)
    if kx == 0.0 && kz == 0.0
        return polarization_ratio(v, 0.0, 1.0)
    else
        return polarization_ratio(v, kx, kz)
    end
end

"""
    handedness(v; kx=0.0, kz=0.0, threshold=1e-5)
    handedness(v, ω; kx=0.0, kz=0.0, threshold=1e-5)
    handedness(v, ω, kx, kz; threshold=1e-5)

Determine the physical polarization handedness of the wave mode represented by the eigenvector `v`
with complex frequency `ω` and wavevector `(kx, kz)` relative to the background magnetic field direction (+z direction).

We assume the wave phase convention is e^{i(k_z z - ω t)} with wave propagation along +z.
This physical handedness correctly projects the electric field perpendicular to the wavevector `k`
to support oblique modes, and factors in the sign of real(ω) to determine the real-time 
rotation sense of the physical field.

Returns:
* `:R` (Right-handed/whistler-like/electron-sense)
* `:L` (Left-handed/Alfvén-like/ion-sense)
* `:linear` (Linearly polarized or near-linear within the threshold)
"""
function handedness(v, ω, kx, kz; threshold=1e-5)
    Ex, Ey, Ez = electric_field(v)
    
    # Calculate transverse field component perpendicular to k in x-z plane
    k = sqrt(kx^2 + kz^2)
    if k == 0.0
        E_perp = Ex
    else
        cos_θ = kz / k
        sin_θ = kx / k
        E_perp = Ex * cos_θ - Ez * sin_θ
    end
    
    abs_E_perp = abs(E_perp)
    abs_Ey = abs(Ey)
    if abs_E_perp < 1e-10 && abs_Ey < 1e-10
        return :linear
    end
    P = Ey / E_perp
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

function handedness(v; kx=0.0, kz=0.0, threshold=1e-5)
    if kx == 0.0 && kz == 0.0
        return handedness(v, 1.0, 0.0, 1.0; threshold)
    else
        return handedness(v, 1.0, kx, kz; threshold)
    end
end

function handedness(v, ω; kx=0.0, kz=0.0, threshold=1e-5)
    if kx == 0.0 && kz == 0.0
        return handedness(v, ω, 0.0, 1.0; threshold)
    else
        return handedness(v, ω, kx, kz; threshold)
    end
end

