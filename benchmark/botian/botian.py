"""Dump BO-T's per-wavevector eigenvalues and timings to results/.

Only the eigensolve of BO-T's reduced matrix ("dynamic roots") is run: the structural
poles it appends afterwards and its root selector have no PlasmaBO counterpart.
"""

import sys
from pathlib import Path
from dataclasses import replace
from time import perf_counter

import numpy as np
from scipy.io import loadmat
from threadpoolctl import threadpool_limits
from boarbitrary_py.unified_runner import initialize_registered_case
from boarbitrary_py.legacy_pdrk.scan import iter_legacy_scan_points
from boarbitrary_py.legacy_pdrk.minimal import solve_minimal_point
from boarbitrary_py.scan import iter_scan_points
from boarbitrary_py.hh_coefficients import compute_hh_coefficients
from boarbitrary_py.hh_minimal import solve_hh_minimal

botian = Path(sys.argv[1]).resolve()
out = Path(__file__).parent / "results"
out.mkdir(exist_ok=True)


def legacy(ctx):
    return iter_legacy_scan_points(ctx), lambda p: solve_minimal_point(ctx, p)


def hh(ctx):
    return iter_scan_points(ctx), lambda p: solve_hh_minimal(compute_hh_coefficients(ctx, p, quadrature_order=80))


CASES = {
    "umeda": ("bo2018_ringbeam_umeda12", "umeda_ring_beam", legacy),
    "firehose": ("bo2025_firehose_astfalk17_theta45", "firehose", hh),
}

def dump(name, points, solve):
    solve(points[0])  # warm-up
    with open(out / f"{name}.csv", "w") as f:
        for i, p in enumerate(points):
            t = perf_counter()
            roots = solve(p).dynamic_roots
            t = perf_counter() - t
            for w in roots:
                f.write(f"{i},{p.kx:.17g},{p.kz:.17g},{t:.6g},{w.real:.17g},{w.imag:.17g}\n")
    print(name, "done")


def at(p, ctx, k_norm, theta_deg):
    """Scan point `p` moved to (k_norm, theta_deg)."""
    kz, kx = k_norm * np.cos(np.deg2rad(theta_deg)), k_norm * np.sin(np.deg2rad(theta_deg))
    return replace(p, par=np.array([k_norm, theta_deg, kz, kx, 0, 0]), k_norm=k_norm, theta_deg=theta_deg,
                   kz_norm=kz, kx_norm=kx, k=k_norm * ctx.kn, kz=kz * ctx.kn, kx=kx * ctx.kn)


with threadpool_limits(1, "blas"):
    for name, (case_id, inputs, family) in CASES.items():
        ctx = initialize_registered_case(case_id, workspace=botian / "inputs" / inputs).context
        points, solve = family(ctx)
        points = tuple(points)
        dump(name, points, solve)
        if name == "umeda":  # the (k, θ) instability map
            grid = [at(points[0], ctx, k, th) for k in np.arange(1, 31) / 100 for th in range(5, 90, 5)]
            dump("umeda_map", grid, solve)

# Firehose proton Hermite coefficients, so both codes expand the same distribution.
fvc = loadmat(botian / "inputs/firehose/paper_cases/firehose_astfalk17_theta45/coefficients/fvceff1.mat")["fvc"][0, 0]
np.savetxt(out / "firehose_alm.csv", fvc["alm"], delimiter=",",
           header=f"w_perp={float(fvc['vtp'].item()):.17g} w_para={float(fvc['vtz'].item()):.17g}")
