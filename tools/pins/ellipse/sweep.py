#!/usr/bin/env python
"""DiffractoKit slice 3b, optional check: the clean sweep against SciPy.

Renders the 20 clean sweep rings itself (the planted ring of Tests/DiffractoKitTestSupport/EllipseFixtures.swift with no
noise and no tail, ported line for line to Python floats and stored as Float32; the FNV-1a hash of every ring equals the Swift
renderer's, checked when this file was written), takes fit1D's start the way the test does (the plant's centre + (1.3, -0.9) px,
annulus b - 20 to a + 20), reads the annulus fit1D
reads (radius > inner and <= outer about the start centre, finite pixels), and solves each ring with SciPy 1.18.1's
`least_squares` for the minimum of fit1D's own objective

    sum( ((A dx^2 + B dx dy + C dy^2 - 1) * value)^2 ),   p = (x0, y0, A, B, C), from fit1D's start (1/R0^2, 0, 1/R0^2),

under two parameter scalings that must agree (centre 1e-5 px, ABC 1e-7 relative). The optimum is converted to (a, b, theta)
by NumPy's symmetric eigensolver (py4DSTEM's convention, a along the smaller eigenvalue's axis, theta in [0, pi)) and written,
with the plant, the start, the pixels' FNV-1a hash and the pixel and bin counts, to Tests/DiffractoGeometryTests/Fixtures/scipy-ellipse-sweep-pins.json. The
Swift test (EllipseCalibrationTests, testCleanSweepEqualsSciPysOptimum) fits the same rings with the kit's fit1D and holds it to
the registered bands: centre 1e-5 px, a and b 1e-7 relative, theta 1e-6 rad. This is the port of the sweep of mac4DSTEM
tools/ellipse-calibration-test/reference.py (sweep_optimum, axes_of), on the kit's own rings.

Run with the pinned interpreter:  ../References/python/orix/bin/python -I -B tools/pins/ellipse/sweep.py
"""
import json
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import least_squares

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = ROOT / "Tests" / "DiffractoGeometryTests" / "Fixtures" / "scipy-ellipse-sweep-pins.json"


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


SWEEP_THETAS = [20.0, 45.0, 67.5, 90.0, 112.5, 135.0, 160.0]


def plants():
    """The 20 clean sweep rings: square 112, 128 and 144 px, R 32-52 px, eps 0 / 1 / 3 %, unbroken or a 60 degree arc missing,
    theta from the sweep's list, a deterministic seed each (the seed places the centre and the gap)."""
    out = []
    for size, radius in [(112, 32.0), (128, 40.0), (144, 52.0)]:
        for eps in [0.0, 0.01, 0.03]:
            for gap in [0.0, 60.0]:
                if len(out) < 20:
                    out.append((size, radius, eps, SWEEP_THETAS[len(out) % 7] if eps > 0 else 0.0, gap))
    while len(out) < 20:
        out.append((128, 44.0, 0.03, SWEEP_THETAS[len(out) % 7], 0.0))
    rings = []
    for n, (size, radius, eps, theta, gap) in enumerate(out):
        rings.append({"name": "sweep %d n%d R%d e%s t%s g%d" % (n, size, int(radius), eps, theta, int(gap)), "size": size,
                      "radius": radius, "eps": eps, "thetaDegrees": theta, "gap": gap, "seed": 100 + n})
    return rings


def render(p):
    """EllipsePlant.render with noise 0 and tail 0: I = exp(-((rho - 1) R)^2 / 18), Float32 pixels, row-major (qx the row)."""
    size, radius, eps, seed, gap = p["size"], p["radius"], p["eps"], p["seed"], p["gap"]
    b = 2 * radius / (2 + eps)
    a = b * (1 + eps)
    cx = size / 2 + 0.37 + 2.5 * math.sin(1.3 * seed + 0.2)
    cy = size / 2 - 0.41 + 2.5 * math.cos(1.9 * seed + 0.7)
    t = p["thetaDegrees"] * math.pi / 180
    sin2, cos2 = math.sin(t) * math.sin(t), math.cos(t) * math.cos(t)
    A = sin2 / (b * b) + cos2 / (a * a)
    C = cos2 / (b * b) + sin2 / (a * a)
    B = 2 * (b * b - a * a) * math.sin(t) * math.cos(t) / (a * a * b * b)
    mean = (a + b) / 2
    gap_start = float((47 * seed) % 360)
    pixels = np.zeros(size * size, dtype="<f4")
    for qx in range(size):
        for qy in range(size):
            dx, dy = qx - cx, qy - cy
            rho = math.sqrt(max(0.0, A * dx * dx + B * dx * dy + C * dy * dy))
            distance = (rho - 1) * mean
            value = math.exp(-(distance * distance) / (2 * 3.0 * 3.0))
            if gap > 0:
                angle = math.fmod(math.atan2(dy, dx) * 180 / math.pi - gap_start, 360.0)
                if angle < 0:
                    angle += 360
                if angle < gap:
                    value = 0.0
            pixels[qx * size + qy] = value
    return pixels, cx, cy, a, b


def main():
    rings, spreads = [], []
    for p in plants():
        size = p["size"]
        pixels, cx, cy, a_plant, b_plant = render(p)
        sx, sy, inner, outer = cx + 1.3, cy - 0.9, b_plant - 20, a_plant + 20
        image = pixels.reshape(size, size)
        qx, qy = np.mgrid[0:size, 0:size].astype(np.float64)
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
            "name": p["name"], "size": size, "radius": p["radius"], "eps": p["eps"], "thetaDegrees": p["thetaDegrees"],
            "gap": p["gap"], "seed": p["seed"], "startQX": sx, "startQY": sy, "innerRadius": inner, "outerRadius": outer,
            "pixelHash": fnv(pixels), "sampleCount": int(v.size), "occupiedBins": bins,
            "optimumCenterQX": float(sol[0]), "optimumCenterQY": float(sol[1]),
            "optimumA": a, "optimumB": b, "optimumTheta": theta,
            "scipySpreadABC": spread_abc, "scipySpreadCentre": spread_centre,
        })
    json.dump({"scipy": __import__("scipy").__version__, "numpy": np.__version__, "rings": rings}, open(OUT, "w"), indent=1)
    print("sweep.py: %d rings -> %s (scipy scaling spread at most %.2g relative in ABC, %.2g px in the centre)"
          % (len(rings), OUT, max(s[0] for s in spreads), max(s[1] for s in spreads)))


main()
