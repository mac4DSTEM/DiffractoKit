#!/usr/bin/env python
"""DiffractoKit slice 3b, optional gate (registration docs/registration-slice-3.md section 3): the clean sweep against SciPy.

Reads the 20 clean sweep rings tools/slice-3b/apps-parity/run.sh dumped (the raw Float32 pixels the Swift renderer made,
.build/apps-parity-3b/sweep/sweep-<n>.f32, and their parameters in .build/apps-parity-3b/dump-O.json), reads the annulus fit1D
reads (radius > inner and <= outer about the start centre, finite pixels), and solves each ring with SciPy 1.18.1's
`least_squares` for the minimum of fit1D's own objective

    sum( ((A dx^2 + B dx dy + C dy^2 - 1) * value)^2 ),   p = (x0, y0, A, B, C), from fit1D's start (1/R0^2, 0, 1/R0^2),

under two parameter scalings that must agree (centre 1e-5 px, ABC 1e-7 relative). The optimum is converted to (a, b, theta)
by NumPy's symmetric eigensolver (py4DSTEM's convention, a along the smaller eigenvalue's axis, theta in [0, pi)) and written,
with the pixels' FNV-1a hash and the pixel and bin counts, to Tests/DiffractoGeometryTests/Fixtures/slice-3b-sweep.json. The
Swift test (EllipseCalibrationTests, testCleanSweepEqualsSciPysOptimum) fits the same rings with the kit's fit1D and holds it to
the registered bands: centre 1e-5 px, a and b 1e-7 relative, theta 1e-6 rad. This is the port of the sweep of mac4DSTEM
tools/ellipse-calibration-test/reference.py (sweep_optimum, axes_of), on the kit's own rings.

Run with the pinned interpreter:  ../References/python/orix/bin/python -I -B tools/slice-3b/sweep.py
"""
import json
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import least_squares

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORK = ROOT / ".build" / "apps-parity-3b"
OUT = ROOT / "Tests" / "DiffractoGeometryTests" / "Fixtures" / "slice-3b-sweep.json"


def fnv(pixels):
    h = 0xCBF29CE484222325
    for p in pixels.view("<u4"):
        v = int(p)
        for _ in range(4):
            h = ((h ^ (v & 0xFF)) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
            v >>= 8
    return "%016x" % h


def axes_of(A, B, C):
    """(a, b, theta) from the quadratic form: a along the eigenvector of the smaller eigenvalue, in (qx, qy), mod pi."""
    values, vectors = np.linalg.eigh(np.array([[A, B / 2.0], [B / 2.0, C]]))
    major = vectors[:, 0]
    return (1.0 / math.sqrt(values[0]), 1.0 / math.sqrt(values[1]), math.atan2(major[1], major[0]) % math.pi)


def optimum(x, y, v, sx, sy, inner, outer):
    def residual(p):
        dx, dy = x - p[0], y - p[1]
        return (p[2] * dx * dx + p[3] * dx * dy + p[4] * dy * dy - 1.0) * v

    r = (inner + outer) / 2.0
    p0 = [sx, sy, 1 / r**2, 0.0, 1 / r**2]
    solutions = [least_squares(residual, p0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, x_scale=scale,
                               max_nfev=5000).x
                 for scale in ("jac", [1.0, 1.0, 1e-5, 1e-5, 1e-5])]
    spread_abc = float(np.max(np.abs(solutions[0][2:] - solutions[1][2:])) / np.linalg.norm(solutions[0][2:]))
    spread_centre = float(np.max(np.abs(solutions[0][:2] - solutions[1][:2])))
    if spread_abc > 1e-7 or spread_centre > 1e-5:
        raise SystemExit("sweep.py: scipy's optimum is not stable under rescaling (ABC %g, centre %g px)" % (spread_abc, spread_centre))
    return solutions[0], spread_abc, spread_centre


def main():
    dump = json.load(open(WORK / "dump-O.json"))
    rings, spreads = [], []
    for n, c in enumerate(dump):
        p = c["plant"]
        if p["kind"] != "sweep":
            continue
        size = p["size"]
        pixels = np.fromfile(WORK / "sweep" / ("sweep-%d.f32" % n), dtype="<f4")
        assert pixels.size == size * size
        if fnv(pixels) != c["pixelHash"]:
            raise SystemExit("sweep.py: %s: the f32 file does not hash to the dump's pixelHash" % p["name"])
        image = pixels.reshape(size, size)
        qx, qy = np.mgrid[0:size, 0:size].astype(np.float64)
        sx, sy, inner, outer = c["startQX"], c["startQY"], c["innerRadius"], c["outerRadius"]
        radius = np.hypot(qx - sx, qy - sy)
        inside = (radius > inner) & (radius <= outer) & np.isfinite(image)
        x, y, v = qx[inside], qy[inside], image[inside].astype(np.float64)
        strong = v >= v.min() + 0.2 * (v.max() - v.min())
        angle = np.arctan2(y[strong] - sy, x[strong] - sx) % (2 * math.pi)
        bins = len(np.unique(np.minimum(35, (angle / (2 * math.pi) * 36).astype(int))))
        if bins < 30:
            raise SystemExit("sweep.py: %s holds %d of 36 bins, below the 30 fit1D needs" % (p["name"], bins))
        sol, spread_abc, spread_centre = optimum(x, y, v, sx, sy, inner, outer)
        spreads.append((spread_abc, spread_centre))
        a, b, theta = axes_of(*sol[2:])
        rings.append({
            "name": p["name"], "index": n, "pixelHash": c["pixelHash"], "sampleCount": int(v.size), "occupiedBins": bins,
            "optimumCenterQX": float(sol[0]), "optimumCenterQY": float(sol[1]),
            "optimumA": a, "optimumB": b, "optimumTheta": theta,
            "scipySpreadABC": spread_abc, "scipySpreadCentre": spread_centre,
        })
    json.dump({"scipy": __import__("scipy").__version__, "numpy": np.__version__, "rings": rings}, open(OUT, "w"), indent=1)
    print("sweep.py: %d rings -> %s (scipy scaling spread at most %.2g relative in ABC, %.2g px in the centre)"
          % (len(rings), OUT, max(s[0] for s in spreads), max(s[1] for s in spreads)))


main()
