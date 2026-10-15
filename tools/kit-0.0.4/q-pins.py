"""DiffractoKit 0.0.4, G2 (docs/registration-0.0.4.md section 3): the pin generator for the exact q from a ring radius and its inverse.

Run, from any directory, with the pinned oracle:

    References/pymatgen-venv/bin/python -I -B tools/kit-0.0.4/q-pins.py [out.json]

It writes Tests/DiffractoGeometryTests/Fixtures/q-pins.json BEFORE the Swift tests read it. Deterministic: no clock, no random
numbers; two runs give the same bytes. Reads slice 3a's Fixtures/slice-3-ring-pins.json (320 cases) for the 1/d identity.

The two closed forms (q in 1/Angstrom, no 2 pi; R and L in pixels; lambda in Angstrom):
  q(R; L, lambda) = 2 sin(atan(R/L) / 2) / lambda                 the Bragg angle theta = atan(R/L)/2, q = 2 sin(theta)/lambda = 1/d
  R(q; L, lambda) = L tan(2 asin(q lambda / 2))                   its inverse, needs q lambda / 2 < 1/sqrt(2) (2 theta < 90 degrees)

TRUTH. Python `decimal` at 60 digits from the exact float64 inputs (lambda is the float64 nearest to the CODATA-2018 relativistic
formula in decimal, exactly as slice 3a's tools/slice-3/ring-pins.py computes it: `dec_wavelength_a`, imported from that file's
formulas, copied here verbatim), then rounded ONCE to float64. Two independent decimal evaluations, which must agree to 1e-50
relative before anything is written:
  (a) ALGEBRAIC, no trigonometric library and no series: with x = R/L and c = cos(atan x) = 1/sqrt(1 + x^2),
      sin(atan(x)/2) = sqrt((1 - c)/2);  with s = q lambda/2 = sin(phi), tan(2 phi) = 2 s sqrt(1 - s^2) / (1 - 2 s^2).
  (b) SERIES: atan by argument halving (atan x = 2 atan(x / (1 + sqrt(1 + x^2))), five times) and its Taylor series, sin and cos by
      their Taylor series, asin as atan(s / sqrt(1 - s^2)), tan as sin / cos.
NumPy 2.5.4's float64 route (np.sin(np.arctan(R / L) / 2) * 2 / lambda; L * np.tan(2 * np.arcsin(q lambda / 2))) is stored beside
the truth (`numpy_*`), and the script fails if it leaves 1e-12 of the truth.

THE GRID (as the registration): V in {60, 80, 120, 200, 300} kV, L_px in {500, 2000, 14082.307097295722, 1e5}, R in {1, 40, 150,
600, 2000} px (100 forward cases); the inverse on the q those give (100) plus q at 0.99 of the Bragg limit, q lambda / 2 =
0.99 / sqrt(2) (20: every V and L); THE IDENTITY: slice 3a's 320 cases, q(R, L_px(R, 2 theta), lambda) = 1/d (stored with the
decimal 1/d and the decimal q of the same float64 inputs, which must agree to 1e-12: L_px is rounded to float64 there);
THE NaN PLANTS: zero, negative, NaN, infinite R / L / lambda / q; q lambda / 2 at the limit (the float64 nearest 1/sqrt(2), as
`maximumLambdaOver2d`), one ulp above it, 0.8, 1, 10 (beyond); an overflow (q = 2 sin(...)/lambda past Double.greatestFiniteMagnitude with
a subnormal lambda; R = L tan(...) past it with L = 1e308) and an underflow (a positive radius so small that atan(R/L) is 0, a
positive q so small that the radius is 0); and one finite plant just inside the limit (one ulp below it), which must be a
finite positive number (its value is ill conditioned by 1/(2e-16) and is not compared to the truth).

The fixture holds {"numpy", "tolerance", "forward": [{name, kv, volts, lambda_A, l_px, r_px, q, numpy_q}], "inverse": [{name, kv,
volts, lambda_A, l_px, q, r_px, numpy_r, kind}], "identity": [{name, d_A, r_px, l_px, lambda_A, inv_d, q}], "nan": [{name, fn,
a, l_px, lambda_A, note}] }. A number that is not finite is written as the string "nan" / "inf" / "-inf".
Exit 0: every self-check holds. 1: one failed (printed).
"""
import resource
# Suite rule (owner, 2026-10-15 12:30): a hard 8 GB memory cap before anything else.
resource.setrlimit(resource.RLIMIT_AS, (8 << 30, 8 << 30))
resource.setrlimit(resource.RLIMIT_DATA, (8 << 30, 8 << 30))
import json
import math
import os
import sys
from decimal import Decimal, getcontext

import numpy as np

NUMPY = "2.5.4"
TOL = 1e-12
getcontext().prec = 60

H_S, C_S, E_S, M0_S = "6.62607015E-34", "299792458", "1.602176634E-19", "9.1093837015E-31"
D_H, D_C, D_E, D_M0 = (Decimal(s) for s in (H_S, C_S, E_S, M0_S))
M_TO_A = Decimal(10) ** 10

VOLTAGES_KV = [60, 80, 120, 200, 300]
CAMERA_LENGTHS_PX = [500.0, 2000.0, 14082.307097295722, 1e5]
RADII_PX = [1.0, 40.0, 150.0, 600.0, 2000.0]
BRAGG_FRACTION = Decimal("0.99")
LIMIT_FLOAT = math.sqrt(0.5)       # Swift: RingCalibration.maximumLambdaOver2d = 0.5.squareRoot()
NAN, INF = float("nan"), float("inf")


def fail(msg):
    print("q-pins.py: FAIL: " + msg, file=sys.stderr)
    sys.exit(1)


def rel(a, b):
    return abs(a - b) / abs(b)


def drel(a, b):
    return abs(a - b) / abs(b)      # Decimals


# ---- decimal closed forms ------------------------------------------------------------------------------------------

def dec_wavelength_a(volts):
    """Copied from tools/slice-3/ring-pins.py: hc / sqrt(eV (2 m0 c^2 + eV)), Angstrom, V in volts."""
    ev = D_E * Decimal(volts)
    return D_H * D_C / (ev * (2 * D_M0 * D_C * D_C + ev)).sqrt() * M_TO_A


def dec_q_algebraic(r, l, lam):
    """(a): q = 2 sqrt((1 - c)/2) / lambda, c = 1 / sqrt(1 + (R/L)^2). No trigonometry."""
    x = Decimal(r) / Decimal(l)
    c = 1 / (1 + x * x).sqrt()
    return 2 * ((1 - c) / 2).sqrt() / Decimal(lam)


def dec_radius_algebraic(q, l, lam):
    """(a): R = L 2 s sqrt(1 - s^2) / (1 - 2 s^2), s = q lambda / 2."""
    s = Decimal(q) * Decimal(lam) / 2
    if not (0 < s * s * 2 < 1):
        fail("dec_radius_algebraic outside 0 < s < 1/sqrt(2): s = %s" % s)
    return Decimal(l) * 2 * s * (1 - s * s).sqrt() / (1 - 2 * s * s)


def _atan_series(x):
    total, term, x2, sign = Decimal(0), x, x * x, 1
    for n in range(0, 80):
        total += sign * term / (2 * n + 1)
        term *= x2
        sign = -sign
    return total


def dec_atan(x):
    """atan by five argument halvings (|x| <= 4 -> ~0.01) then the Taylor series (80 terms of 1e-4^n: far beyond 60 digits)."""
    k = 5
    for _ in range(k):
        x = x / (1 + (1 + x * x).sqrt())
    return _atan_series(x) * (2 ** k)


def dec_sin(x):
    total, term, x2 = Decimal(0), x, x * x
    for n in range(1, 120):
        total += term
        term = -term * x2 / ((2 * n) * (2 * n + 1))
    return total


def dec_cos(x):
    total, term, x2 = Decimal(0), Decimal(1), x * x
    for n in range(1, 120):
        total += term
        term = -term * x2 / ((2 * n - 1) * (2 * n))
    return total


def dec_q_series(r, l, lam):
    """(b): q = 2 sin(atan(R/L)/2)/lambda with the series above."""
    return 2 * dec_sin(dec_atan(Decimal(r) / Decimal(l)) / 2) / Decimal(lam)


def dec_radius_series(q, l, lam):
    """(b): R = L tan(2 asin(s)), asin s = atan(s / sqrt(1 - s^2)), tan = sin/cos."""
    s = Decimal(q) * Decimal(lam) / 2
    phi = dec_atan(s / (1 - s * s).sqrt())
    return Decimal(l) * dec_sin(2 * phi) / dec_cos(2 * phi)


def truth_q(r, l, lam):
    a, b = dec_q_algebraic(r, l, lam), dec_q_series(r, l, lam)
    if drel(a, b) > Decimal("1e-50"):
        fail("q(R=%r, L=%r, lambda=%r): algebraic and series decimal differ by %s" % (r, l, lam, drel(a, b)))
    return a


def truth_radius(q, l, lam):
    a, b = dec_radius_algebraic(q, l, lam), dec_radius_series(q, l, lam)
    if drel(a, b) > Decimal("1e-50"):
        fail("R(q=%r, L=%r, lambda=%r): algebraic and series decimal differ by %s" % (q, l, lam, drel(a, b)))
    return a


# ---- numpy float64 routes ------------------------------------------------------------------------------------------

def np_q(r, l, lam):
    return float(2 * np.sin(np.arctan(np.float64(r) / np.float64(l)) / 2) / lam)


def np_radius(q, l, lam):
    return float(l * np.tan(2 * np.arcsin(np.float64(q) * lam / 2)))


# ---- the fixture ---------------------------------------------------------------------------------------------------

def jsonable(x):
    if isinstance(x, float) and not math.isfinite(x):
        return "nan" if math.isnan(x) else ("inf" if x > 0 else "-inf")
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if isinstance(x, list):
        return [jsonable(v) for v in x]
    return x


HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "..", "..", "Tests", "DiffractoGeometryTests", "Fixtures")
DEFAULT_OUT = os.path.join(FIXTURES, "q-pins.json")
RING_PINS = os.path.join(FIXTURES, "slice-3-ring-pins.json")


def nan_plants():
    """Every plant must be NaN: fn 'q' takes (a = R, l_px, lambda_A); fn 'radius' takes (a = q, l_px, lambda_A)."""
    plants = []

    def add(name, fn, a, l_px, lam, note):
        plants.append({"name": name, "fn": fn, "a": a, "l_px": l_px, "lambda_A": lam, "note": note})

    L, LAM, R, Q = 2000.0, 0.025079340450548003, 150.0, 5.0     # a clean base: 200 kV
    for fn, base in (("q", R), ("radius", 1.0)):
        for tag, bad in (("zero", 0.0), ("negative-zero", -0.0), ("negative", -1.0), ("nan", NAN), ("inf", INF), ("minus-inf", -INF)):
            add("%s-first-%s" % (fn, tag), fn, bad, L, LAM, "the radius / q input is not finite and positive")
            add("%s-length-%s" % (fn, tag), fn, base, bad, LAM, "the camera length is not finite and positive")
            add("%s-lambda-%s" % (fn, tag), fn, base, L, bad, "the wavelength is not finite and positive")
    # the Bragg limit of the inverse: s = q lambda / 2 must be < 1/sqrt(2) (the float64 `maximumLambdaOver2d`)
    # lambda = 2 makes s = q exactly; lambda = 1 and q = 2 s makes s = q / 2 exactly.
    above = float(np.nextafter(LIMIT_FLOAT, 1.0))
    for tag, s in (("at-limit", LIMIT_FLOAT), ("one-ulp-above", above), ("0.8", 0.8), ("1", 1.0), ("10", 10.0), ("1e300", 1e300)):
        add("radius-bragg-%s" % tag, "radius", s, L, 2.0, "q lambda / 2 = %r is not below 1/sqrt(2)" % s)
    add("radius-bragg-at-limit-lambda1", "radius", 2 * LIMIT_FLOAT, L, 1.0, "q lambda / 2 = the limit again, from lambda = 1")
    # overflow and underflow: not representable, so NaN and never infinity or a plausible 0
    add("q-overflow-subnormal-lambda", "q", 1.0, 1.0, 5e-324, "q = 2 sin(pi/8)/5e-324 is past Double.greatestFiniteMagnitude")
    add("q-overflow-lambda-1e-310", "q", 1e300, 1e300, 1e-310, "q ~ 1.5e310")
    add("radius-overflow-L-1e308", "radius", 0.99 * 2 * 0.5 ** 0.5 / 2.0, 1e308, 2.0, "R = 1e308 tan(...) with s at 0.99 of the limit is ~5e309")
    add("radius-overflow-L-max", "radius", 0.5, 1.7976931348623157e308, 2.0, "R = L tan(2 asin(0.5)) = 1.73 L is past the largest double")
    add("q-underflow-R-subnormal", "q", 5e-324, 1e5, 0.025079340450548003, "atan(R/L) is 0 in float64: q would read 0 for a positive radius")
    add("q-underflow-R-over-L", "q", 1e-300, 1e300, 1.0, "R/L = 1e-600 underflows to 0")
    add("radius-underflow-q-subnormal", "radius", 5e-324, 1e5, 0.025079340450548003, "the radius would read 0 for a positive q")
    add("radius-underflow-tiny-L", "radius", 1e-5, 5e-324, 0.025, "L tan(...) < the smallest subnormal")
    return plants


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    if np.__version__ != NUMPY:
        fail("numpy %s, the lock pins %s" % (np.__version__, NUMPY))

    worst = {"np_q": 0.0, "np_radius": 0.0, "np_roundtrip": 0.0, "np_identity": 0.0, "dec_identity": 0.0}
    forward, inverse = [], []
    for kv in VOLTAGES_KV:
        volts = kv * 1000.0
        lam = float(dec_wavelength_a(volts))
        for l in CAMERA_LENGTHS_PX:
            for r in RADII_PX:
                t = truth_q(r, l, lam)
                q = float(t)
                npq = np_q(r, l, lam)
                worst["np_q"] = max(worst["np_q"], rel(npq, q))
                forward.append({"name": "q-%dkV-L%g-R%g" % (kv, l, r), "kv": kv, "volts": volts, "lambda_A": lam, "l_px": l,
                                "r_px": r, "q": q, "numpy_q": npq})
                # the inverse on this q (a float64 input): truth from the exact float64 q, not from the decimal t
                tr = truth_radius(q, l, lam)
                rr = float(tr)
                npr = np_radius(q, l, lam)
                worst["np_radius"] = max(worst["np_radius"], rel(npr, rr))
                worst["np_roundtrip"] = max(worst["np_roundtrip"], rel(npr, r))
                if drel(tr, Decimal(r)) > Decimal("1e-15"):
                    fail("round trip in decimal: R(q(R)) = %s for R = %r (q rounded to float64)" % (tr, r))
                inverse.append({"name": "radius-%dkV-L%g-R%g" % (kv, l, r), "kind": "roundtrip", "kv": kv, "volts": volts,
                                "lambda_A": lam, "l_px": l, "q": q, "r_px": rr, "numpy_r": npr, "start_r_px": r})
            # q at 0.99 of the Bragg limit: s = q lambda / 2 = 0.99 / sqrt(2)
            q99 = float(BRAGG_FRACTION * Decimal("0.5").sqrt() * 2 / Decimal(lam))
            tr = truth_radius(q99, l, lam)
            npr = np_radius(q99, l, lam)
            s99 = Decimal(q99) * Decimal(lam) / 2
            if abs(float(s99 * Decimal(2).sqrt()) - 0.99) > 1e-15:
                fail("q at 0.99 of the limit: s sqrt 2 = %s" % (s99 * Decimal(2).sqrt()))
            worst["np_radius"] = max(worst["np_radius"], rel(npr, float(tr)))
            inverse.append({"name": "radius-%dkV-L%g-bragg0.99" % (kv, l), "kind": "bragg99", "kv": kv, "volts": volts,
                            "lambda_A": lam, "l_px": l, "q": q99, "r_px": float(tr), "numpy_r": npr, "start_r_px": None})
    if len(forward) != 100 or len(inverse) != 120:
        fail("%d forward, %d inverse cases" % (len(forward), len(inverse)))

    # ---- the identity on slice 3a's 320 cases ----
    with open(RING_PINS) as f:
        ring = json.load(f)
    if len(ring["cases"]) != 320:
        fail("slice-3-ring-pins.json holds %d cases" % len(ring["cases"]))
    identity = []
    for c in ring["cases"]:
        r, l, lam, d = c["r_px"], c["l_px"], c["lambda_A"], c["d_A"]
        t = truth_q(r, l, lam)
        inv_d = Decimal(1) / Decimal(d)
        worst["dec_identity"] = max(worst["dec_identity"], float(drel(t, inv_d)))
        worst["np_identity"] = max(worst["np_identity"], rel(np_q(r, l, lam), float(inv_d)))
        identity.append({"name": c["name"], "d_A": d, "r_px": r, "l_px": l, "lambda_A": lam, "inv_d": float(inv_d), "q": float(t)})
    if worst["dec_identity"] > TOL:
        fail("the identity q(R, L_px, lambda) = 1/d fails in decimal: %.3e" % worst["dec_identity"])

    plants = nan_plants()
    if len({p["name"] for p in plants}) != len(plants):
        fail("duplicate plant names")
    for k in ("np_q", "np_radius", "np_roundtrip", "np_identity"):
        if worst[k] > TOL:
            fail("numpy float64 route leaves the band: %s %.3e" % (k, worst[k]))

    # the plants in the domain's interior must really have a finite truth, and the NaN ones none: spot-check the inputs
    above = float(np.nextafter(LIMIT_FLOAT, 1.0))
    below = float(np.nextafter(LIMIT_FLOAT, 0.0))
    if not (Decimal(below) < Decimal("0.5").sqrt() < Decimal(LIMIT_FLOAT) < Decimal(above)):
        fail("the float64 neighbours of 1/sqrt(2) are not where the Bragg plants assume")
    just_inside = {"name": "radius-just-inside-limit", "q": below, "l_px": 2000.0, "lambda_A": 2.0}

    doc = {"numpy": NUMPY, "tolerance": TOL, "limit_float": LIMIT_FLOAT, "forward": forward, "inverse": inverse,
           "identity": identity, "nan": plants, "just_inside": just_inside}
    with open(out, "w") as f:
        json.dump(jsonable(doc), f, indent=1, sort_keys=True, allow_nan=False)
        f.write("\n")
    print("q-pins: %d forward, %d inverse (%d round trip + %d at 0.99 of the limit), %d identity, %d NaN plants"
          % (len(forward), len(inverse), 100, 20, len(identity), len(plants)))
    print("q-pins: numpy float64 against the decimal truth, worst relative: " + ", ".join("%s %.2e" % (k, v) for k, v in worst.items()))


if __name__ == "__main__":
    main()
