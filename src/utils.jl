"""
    solve_with_threads(f, nthreads)

Execute function `f` with specified BLAS thread count, restoring previous setting afterward.
"""
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
