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
    return @showprogress dt = dt desc = desc for id in carts
        ik, iθ = Tuple(id)
        k = ks[ik]
        θ = θs[iθ]
        kx, kz = k .* sincos(θ)
        f(ik, iθ, kx, kz)
    end
end
