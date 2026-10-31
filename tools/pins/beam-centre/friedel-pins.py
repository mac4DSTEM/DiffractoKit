"""Pins for FriedelCentre (DiffractoGeometry): py4DSTEM f050d207 `get_origin_friedel` (process/calibration/origin.py,
line 332) on single synthetic patterns whose central beam is blocked, with the truth they were planted with.

Run (the py4DSTEM venv, see tools/pins/requirements/py4dstem.txt):
    ../References/python/py4dstem/bin/python -I -B tools/pins/beam-centre/friedel-pins.py [out.json]
Writes Tests/DiffractoGeometryTests/Fixtures/py4dstem-friedel-centre-pins.json (deterministic: seeded generator, no clock).

Each case stores the Float32 pixels as written (the Swift test reads them, it does not regenerate them), the mask when
there is one (true in the pattern, false under the beam stop), the planted centre and py4DSTEM's (qx0, qy0) = (row,
column) with upsample_factor 1 (the half-pixel parabola). py4DSTEM states no tolerance for the function and has no
test of it, so the Swift side holds assert_allclose defaults for float64 (rtol 1e-7).

The diagnostic block prints what py4DSTEM's simple finder `get_origin_single_dp` (the brightest smoothed pixel and the
centre of mass about it) returns on the same patterns: on a blocked beam it lands on a Bragg spot or the beam stop edge.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from py4DSTEM import DataCube
from py4DSTEM.process.calibration.origin import get_origin_friedel, get_origin_single_dp

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Tests/DiffractoGeometryTests/Fixtures/py4dstem-friedel-centre-pins.json"


def pattern(rows, cols, centre, seed, stop_radius, arm, saturate=None, noise=0.01):
    """A spot pattern (two lattice vectors, Friedel symmetric about `centre`), a diffuse halo and a ring, the central beam
    blocked by a disc of `stop_radius` with an arm of half width 1.5 px from the centre to the first-row edge (`arm`), or
    clipped at `saturate`. Multiplicative uniform noise of `noise`."""
    cr, cc = centre
    r, c = np.meshgrid(np.arange(rows, dtype=np.float64), np.arange(cols, dtype=np.float64), indexing="ij")
    dr, dc = r - cr, c - cc
    rho = np.hypot(dr, dc)
    img = 20.0 + 400.0 * np.exp(-rho / 9.0)
    ring = 0.36 * min(rows, cols)
    img += 150.0 * np.exp(-((rho - ring) ** 2) / (2 * 1.3 ** 2))
    g1, g2 = np.array([7.3, 2.1]), np.array([-1.7, 8.4])
    for h in range(-4, 5):
        for k in range(-4, 5):
            if h == 0 and k == 0:
                continue
            g = h * g1 + k * g2
            amp = 3000.0 * np.exp(-(g @ g) / (2 * 14.0 ** 2))
            img += amp * np.exp(-((dr - g[0]) ** 2 + (dc - g[1]) ** 2) / (2 * 1.1 ** 2))
    img += 60000.0 * np.exp(-(rho ** 2) / (2 * 2.2 ** 2))
    rng = np.random.default_rng(seed)
    img *= 1 + noise * rng.uniform(-1, 1, img.shape)
    mask = np.ones((rows, cols), dtype=bool)
    if stop_radius > 0:
        mask &= rho > stop_radius
        if arm:
            mask &= ~((np.abs(dc) < 1.5) & (dr < 0))
        img = np.where(mask, img, 0.0)
    if saturate is not None:
        img = np.minimum(img, saturate)
    return img.astype(np.float32), mask


CASES = [
    # name, rows, cols, centre (row, col), seed, stop radius, arm, saturate, use mask
    ("stop-arm-nomask", 64, 64, (31.3, 33.6), 1, 6.0, True, None, False),
    ("stop-arm-mask", 64, 64, (31.3, 33.6), 1, 6.0, True, None, True),
    ("stop-disc-nomask-rect", 56, 72, (27.8, 36.45), 2, 5.0, False, None, False),
    ("stop-disc-mask-rect", 56, 72, (27.8, 36.45), 2, 5.0, False, None, True),
    ("saturated-centre", 64, 64, (32.5, 30.2), 3, 0.0, False, 900.0, False),
    ("open-beam", 64, 64, (30.75, 32.25), 4, 0.0, False, None, False),
]


def main():
    cases = []
    for name, rows, cols, centre, seed, stop, arm, sat, use_mask in CASES:
        img, mask = pattern(rows, cols, centre, seed, stop, arm, sat)
        nonzero = img[img != 0]
        assert nonzero.min() > 16 and float(img.astype(np.float64).sum()) < 2.0 ** 24, name   # the exact-sum claim below
        dc = DataCube(data=img[None, None, :, :].astype(np.float64))
        qx0, qy0 = get_origin_friedel(dc, mask=mask if use_mask else None, upsample_factor=1)
        qx, qy = float(qx0[0, 0]), float(qy0[0, 0])
        simple = get_origin_single_dp(img.astype(np.float64), 2.2)
        print(f"{name:24s} truth ({centre[0]:7.3f}, {centre[1]:7.3f})  friedel ({qx:8.4f}, {qy:8.4f})  "
              f"single_dp ({float(simple[0]):8.3f}, {float(simple[1]):8.3f})")
        case = {
            "name": name, "rows": rows, "columns": cols,
            "truthRow": centre[0], "truthColumn": centre[1],
            # the shortest decimal that reads back to the same Float32. pixelSum, the exact sum of the Float32 values
            # (every nonzero pixel is above 16, so its last bit is at least 2^-19, and the sum stays below 2^24: about 43
            # bits, so a float64 sum is exact in any order), lets the Swift test prove it read every pixel exactly
            "pixels": [float(str(v)) for v in img.ravel()],
            "pixelSum": math.fsum(float(v) for v in img.ravel()),
            "py4dstemRow": qx, "py4dstemColumn": qy,
        }
        if use_mask:
            case["mask"] = [int(v) for v in mask.ravel()]     # 1 in the pattern, 0 under the beam stop
        cases.append(case)
    doc = {
        "reference": "py4DSTEM f050d207 process/calibration/origin.py get_origin_friedel(datacube, mask, upsample_factor=1)",
        "tolerance": "rtol 1e-7 (assert_allclose float64 defaults; upstream has no test of the function)",
        "cases": cases,
    }
    OUT.write_text(json.dumps(doc, separators=(",", ":")) + "\n")
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
