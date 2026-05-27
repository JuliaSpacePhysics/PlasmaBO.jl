# Compare per-call cost of integer-order J_n(x) across:
#   __besselj_integer        (generic series for |x|<7, quadgk otherwise — kept
#                             here since it is too slow for production use)
#   Bessels.besselj          (Float32/Float64 fast path)
#   SpecialFunctions.besselj (libopenspecfun; BigFloat / Double64 via ext)
#
# Ratio column = __besselj_integer / reference (Bessels for Float64, SpecFun otherwise)

using BenchmarkTools
using Bessels
using DoubleFloats: Double64
using Logging
using MultiFloats
using Printf
using QuadGK: quadgk
using SpecialFunctions

# MultiFloats installs transcendentals dynamically; silence overwrite warnings.
with_logger(NullLogger()) do
    MultiFloats.use_bigfloat_transcendentals()
end

# Generic fallback: power series for |x|<7, integral representation otherwise.
# Dropped from PlasmaBO because the quadgk branch is unusably slow for extended
# precision; retained here purely as a benchmark baseline.
function __besselj_integer(n::Integer, x::T) where {T <: Real}
    if abs(x) < T(7)
        n < 0 && return isodd(n) ? -__besselj_integer(-n, x) : __besselj_integer(-n, x)
        term = one(T)
        for k in 1:n
            term *= x / (T(2) * T(k))
        end
        value = term
        z = -(x * x) / T(4)
        for k in 1:10_000
            term *= z / (T(k) * T(n + k))
            next = value + term
            (next == value || abs(term) <= eps(one(T)) * max(one(T), abs(next))) && return next
            value = next
        end
        error("Bessel J series failed to converge for n=$n and x=$x")
    end
    value = quadgk(
        θ -> cos(T(n) * θ - x * sin(θ)), zero(T), T(π);
        rtol = sqrt(eps(one(T))),
    )[1]
    return value / T(π)
end

# Probes `f(n, x)` once with stderr silenced (some unsupported combos overflow
# the stack and the runtime prints to stderr before the catch can fire).
# Returns the median benchmark time, or `nothing` if the probe errors.
function bench_ns(f, n, x)
    try
        redirect_stderr(devnull) do
            f(n, x)
        end
    catch
        return nothing
    end
    return median(@benchmark $f($n, $x) samples = 20 evals = 10).time
end

fmt_ns(::Nothing) = @sprintf("%13s", "N/A")
fmt_ns(t::Real) = @sprintf("%10.1f ns", t)

fmt_ratio(::Nothing, _) = @sprintf("%8s", "N/A")
fmt_ratio(_, ::Nothing) = @sprintf("%8s", "N/A")
fmt_ratio(num::Real, den::Real) = @sprintf("%7.2fx", num / den)

const NS = (0, 1, 2, 5, 10)
const XS = (0.1, 1.0, 3.5, 7.0, 15.0)

# (label, x-converter, columns); the last column is the ratio reference.
const SECTIONS = [
    (
        "Float64", identity, [
            ("__besselj", __besselj_integer),
            ("Bessels", Bessels.besselj),
            ("SpecFun", SpecialFunctions.besselj),
        ], 2,
    ),  # ratio denominator = column 2 (Bessels)
]

for (label, conv, cols, ref_idx) in SECTIONS
    println("\n=== $label ===")
    @printf("%-4s %-6s", "n", "x")
    for (h, _) in cols
        @printf(" %12s", h)
    end
    ref_idx > 0 && @printf("   __besselj/%s", cols[ref_idx][1])
    println()
    for n in NS, x in XS
        xv = conv(x)
        @printf("%-4d %-6.2f", n, x)
        ts = [bench_ns(f, n, xv) for (_, f) in cols]
        for t in ts
            print(" ", fmt_ns(t))
        end
        ref_idx > 0 && print("   ", fmt_ratio(ts[1], ts[ref_idx]))
        println()
    end
end

# Side-by-side comparison of the two ~106-bit double-double formats.
# Bessels.besselj / SpecialFunctions.besselj are tried for both; missing
# methods print "N/A" instead of erroring.
println("\n=== Extended precision: Double64 vs Float64x2 (~106-bit) ===")
@printf(
    "%-4s %-6s %13s %13s %13s   %13s %13s %13s   %8s\n",
    "n", "x",
    "D64.__bj", "D64.Bessels", "D64.SpecFun",
    "F64x2.__bj", "F64x2.Bessels", "F64x2.SpecFun",
    "F64x2/D64"
)
for n in NS, x in XS
    xd = Double64(x)
    xm = Float64x2(x)
    td_bj = bench_ns(__besselj_integer, n, xd)
    td_bs = bench_ns(Bessels.besselj, n, xd)
    td_sf = bench_ns(SpecialFunctions.besselj, n, xd)
    tm_bj = bench_ns(__besselj_integer, n, xm)
    tm_bs = bench_ns(Bessels.besselj, n, xm)
    tm_sf = bench_ns(SpecialFunctions.besselj, n, xm)
    @printf("%-4d %-6.2f", n, x)
    for t in (td_bj, td_bs, td_sf, tm_bj, tm_bs, tm_sf)
        print(" ", fmt_ns(t))
    end
    print("   ", fmt_ratio(tm_bj, td_bj))
    println()
end
