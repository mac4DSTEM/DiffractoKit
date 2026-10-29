"""DiffractoKit slice 4: the pin generator for the orientation and phase-matching module.

Run, from any directory, under the memory cap, with the pymatgen venv (Python 3.14, NumPy 2.5.4, SciPy 1.18.1):

    References/python/orix/bin/python -I -B DiffractoKit/tools/slice-4/pins.py

It writes Tests/DiffractoOrientationTests/Fixtures/slice-4-pins.json BEFORE any Swift of slice 4 is trusted. Deterministic: no
clock, seeded random numbers (np.random.default_rng), sorted keys, so two runs give the same bytes.

Two interpreters, one memory cap: this script runs the py4DSTEM half (tools/slice-4/pins-py4dstem.py) as a subprocess with
References/python/py4dstem/bin/python (Python 3.12, py4DSTEM f050d207; read-only, never installed into) and reads its JSON
from stdout. SciPy 1.18.1 and the NumPy oracle are here, because the py4DSTEM venv carries SciPy 1.18.0.

What is an oracle for what
--------------------------
V2  "reference": NumPy float64 written from the DEFINITIONS in this file, the structure factors from the
    py4DSTEM electron-scattering table (slice 2's sha-locked file) with the cell/g conventions py4DSTEM uses, every hkl with
    |F| > 1e-4 (strict, slice 2) and |g| <= kMax enumerated by a bound that cannot miss a lattice point, the FLAT excitation
    s = g.n, I = |F|^2 exp(-(s/w)^2), the in-plane frame (e1, e2) of `ZoneAxisFrame`. This is what the Swift library must equal
    (G4, and the pre-cap entries of G2).
V0  py4DSTEM f050d207 as it is called (Crystal.generate_diffraction_pattern at a huge accelerating voltage, sigma = w/sqrt(2),
    tol_excitation_error_mult = slab/sigma, tol_intensity 0): its hkl set, g in its projected frame, relative intensity.
V1  NumPy flat excitation on py4DSTEM's OWN hkl, g and |F|^2 (float32): V0 against V1 is the lambda -> 0 residual;
    V0 against V2 is py4DSTEM's float32 structure-factor floor plus any hkl py4DSTEM's tiling misses.
V200 py4DSTEM at 200 kV with the curved Ewald sphere: a parity row (hkl entering and leaving the slab), not a gate.

Exit status: 0 all self-checks and all predicted counts hold; 1 a self-check failed (assert, loud); 3 the file was written but a
registered or predicted number differs from the measured one (printed twice, never "fixed" here).
"""
import hashlib
import json
import math
import os
import subprocess
import sys
import warnings
from importlib.metadata import version as dist_version

import numpy as np
import scipy
from scipy.spatial.transform import Rotation

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser

# ----------------------------------------------------------------------------------------------------------------
# Paths (relative to this script, so it runs from any directory)
# ----------------------------------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", ".."))
SUITE = os.path.abspath(os.path.join(KIT, ".."))
TABLE_PATH = os.path.join(KIT, "References", "py4dstem-scattering", "scattering_factors.txt")
PY4D_PYTHON = os.path.join(SUITE, "References", "python", "py4dstem", "bin", "python")
PY4D_HELPER = os.path.join(HERE, "pins-py4dstem.py")
BETA_CIF = os.path.join(KIT, "References", "crystal-structures", "beta_double_prime_Mg5Si6.cif")
SLICE2_PINS = os.path.join(KIT, "Tests", "DiffractoScatteringTests", "Fixtures", "slice-2-pins.json")
OUT_DIR = os.path.join(KIT, "Tests", "DiffractoOrientationTests", "Fixtures")
OUT = os.path.join(OUT_DIR, "slice-4-pins.json")

# ----------------------------------------------------------------------------------------------------------------
# The registered constants (the defaults and the method's constants)
# ----------------------------------------------------------------------------------------------------------------
TABLE_SHA256 = "be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4"
K_MAX = 1.6
SLAB = 0.05
WIDTH = 0.03
FLOOR = 0.05
CAP = 48
IN_PLANE_STEP = 2.0
CUT = 1e-4                  # |F| cut, strict (slice 2)
TIE_TOL = 1e-12             # slice 2 / slice 1 3.8: the kit's reflection order
Q_FLOOR = 1e-9              # the length floor of a projected vector
BASIS_SWITCH = 0.9          # |n_x| < 0.9: reference axis x, else y
SEED = 20261015

CHECKS = []
PREDICTIONS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    assert ok, f"SELF-CHECK FAILED: {name} {detail}"


def predict(name, predicted, measured):
    PREDICTIONS.append((name, predicted, measured))
    if predicted != measured:
        print(f"!!! PREDICTION DIFFERS: {name}: registered {predicted}, measured {measured}", flush=True)


# ----------------------------------------------------------------------------------------------------------------
# f_e, lattice, reflections (slice 2's definitions, copied, not imported)
# ----------------------------------------------------------------------------------------------------------------
def load_table():
    data = open(TABLE_PATH, "rb").read()
    got = hashlib.sha256(data).hexdigest()
    check("the scattering table's sha256 is the lock", got == TABLE_SHA256, got)
    t = np.loadtxt(TABLE_PATH, dtype=np.float64)
    check("the table is 103 x 10", t.shape == (103, 10))
    return t


def f_e(table, z, g2):
    row = table[z - 1]
    a, b = row[0:10:2], row[1:10:2]
    g2 = np.asarray(g2, dtype=np.float64)
    return sum(a[i] * (2 + b[i] * g2) / (1 + b[i] * g2) ** 2 for i in range(5))


def lattice(cell):
    a, b, c, al, be, ga = cell
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    sg = math.sin(math.radians(ga))
    vol_formula = a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    A = np.array([[a, 0, 0],
                  [b * cg, b * sg, 0],
                  [c * cb, -c * (cb * cg - ca) / sg, vol_formula / (a * b * sg)]])
    vol = abs(np.linalg.det(A))
    check(f"cell{cell}: |det A| equals the closed-form volume to 1e-12 relative", abs(vol - vol_formula) <= 1e-12 * vol_formula)
    L = Lattice.from_parameters(a, b, c, al, be, ga, vesta=True)
    check(f"cell{cell}: A equals pymatgen's vesta matrix to 1e-12", np.allclose(A, np.array(L.matrix), rtol=0, atol=1e-12 * max(a, b, c)))
    B = np.linalg.inv(A).T          # a_i . b_j = delta_ij, g = h b1 + k b2 + l b3, no 2 pi
    check("a_i . b_j = delta_ij", np.allclose(A @ B.T, np.eye(3), rtol=0, atol=1e-12))
    return A, B, float(vol)


def order_reflections(H, g):
    """Increasing |g|; equal |g| (within 1e-12 relative, taken from the cluster's first) by descending (h, k, l)."""
    idx = np.argsort(g, kind="stable")
    out = []
    n = len(idx)
    i = 0
    while i < n:
        g0 = g[idx[i]]
        j = i + 1
        while j < n and g[idx[j]] <= g0 * (1 + TIE_TOL):
            j += 1
        out.extend(sorted(idx[i:j].tolist(), key=lambda m: tuple(int(x) for x in H[m]), reverse=True))
        i = j
    return out


def enumerate_box(A, B, kmax):
    """Every (h, k, l) != 0 with |g| <= kmax. h = g . a_1 for a lattice point g, so |h| <= kmax |a_1|: the box
    ceil(kmax |a_i|) + 1 cannot miss one (py4DSTEM's tiling can: see the per-phase record)."""
    n = [int(math.ceil(kmax * np.linalg.norm(A[i]))) + 1 for i in range(3)]
    axes = [np.arange(-m, m + 1) for m in n]
    H = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3).astype(np.int64)
    H = H[np.any(H != 0, axis=1)]
    g = np.linalg.norm(H.astype(float) @ B, axis=1)
    keep = g <= kmax
    return H[keep], g[keep]


def structure_factors(table, vol, H, glen, sites):
    """F(hkl) = (1/V) sum_sites f_e(Z, |g|^2) occ exp(2 pi i (h x + k y + l z)); sites = [(Z, (x, y, z), occ)]."""
    Hf = H.astype(float)
    F = np.zeros(len(H), dtype=np.complex128)
    for (z, pos, occ) in sites:
        F += f_e(table, z, glen ** 2) * occ * np.exp(2j * np.pi * (Hf @ np.array(pos, dtype=float)))
    return F / vol


# ----------------------------------------------------------------------------------------------------------------
# The phases
# ----------------------------------------------------------------------------------------------------------------
FCC = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
HCP = [(0.0, 0.0, 0.0), (2.0 / 3.0, 1.0 / 3.0, 0.5)]
CUBIC_AXES = [[0, 0, 1], [0, 1, 1], [1, 1, 1], [1, 1, 2], [1, 2, 3]]


def wrap_unit(x):
    v = x % 1.0
    return 0.0 if v >= 1.0 - 1e-12 else float(v)


def beta_double_prime_sites():
    """The apps' Crystal.betaDoublePrime (mac4DSTEM 7d6fdbd1, Crystal.swift): Andersen et al. 1998 Table 3 set 3 under C2/m,
    in the apps' expansion order and expressions: Mg 2a x2, then 4i Mg2, Mg3, Si1, Si2, Si3, each (x,0,z), (-x,0,-z),
    (x+1/2,1/2,z), (-x+1/2,1/2,-z). Not wrapped into [0, 1), as the apps."""
    def four_i(z, x, zc):
        return [(z, (x, 0.0, zc), 1.0), (z, (-x, 0.0, -zc), 1.0), (z, (x + 0.5, 0.5, zc), 1.0), (z, (-x + 0.5, 0.5, -zc), 1.0)]
    mg, si = 12, 14
    s = [(mg, (0.0, 0.0, 0.0), 1.0), (mg, (0.5, 0.5, 0.0), 1.0)]
    s += four_i(mg, 0.3459, 0.089) + four_i(mg, 0.430, 0.652)
    s += four_i(si, 0.0565, 0.649) + four_i(si, 0.1885, 0.224) + four_i(si, 0.2171, 0.617)
    return s


def beta_cif_check(sites):
    """The apps' 22 sites equal, modulo 1 and as a set, what pymatgen expands from the beta'' CIF the kit's slice 2 pinned."""
    if not os.path.isfile(BETA_CIF):
        return {"checked": False, "reason": "References/crystal-structures/beta_double_prime_Mg5Si6.cif absent"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = CifParser(BETA_CIF, frac_tolerance=0.0).parse_structures(primitive=False)[0]
    cif = sorted((int(s.specie.Z), tuple(wrap_unit(x) for x in s.frac_coords)) for s in st)
    app = sorted((z, tuple(wrap_unit(x) for x in pos)) for (z, pos, _) in sites)
    check("beta'' CIF expansion has 22 sites", len(cif) == 22 == len(app))
    worst = 0.0
    for (zc, pc), (za, pa) in zip(cif, app):
        # sorted order may differ in the last digits; match by nearest instead
        pass
    unmatched = 0
    for (za, pa) in app:
        best = min((max(min(abs(a - b), 1 - abs(a - b)) for a, b in zip(pa, pc)) for (zc, pc) in cif if zc == za), default=9.0)
        worst = max(worst, best)
        if best > 1e-9:
            unmatched += 1
    check("beta'': every app site (modulo 1) is a CIF site to 1e-9", unmatched == 0, f"worst {worst}")
    return {"checked": True, "file": "beta_double_prime_Mg5Si6.cif", "sha256": hashlib.sha256(open(BETA_CIF, "rb").read()).hexdigest(),
            "worst_site_difference_mod1": worst}


def theta_prime_sites():
    """theta'-Al2Cu: a = 4.04, c = 5.80 A, body-centred tetragonal, Cu 2a (0,0,0) + (1/2,1/2,1/2), Al 4d
    (0,1/2,1/4), (1/2,0,1/4), (0,1/2,3/4), (1/2,0,3/4): the fluorite-derived arrangement of Silcock 1953 (Cu on the Ca sites,
    Al on the F sites of CaF2 seen along its cube diagonal-free tetragonal cell). NO app or kit definition exists (the apps' tests
    carry none, the kit has no Al2Cu CIF): the sites are written out here, not expanded from a space group, because pymatgen's
    I-4m2 (mirrors on the axes) holds the site (0,1/2,1/4) as a 2-fold orbit, so the 4d orbit needs I-42m or the sites listed;
    the pin carries the sites so the Swift builds exactly this structure."""
    sites = [(29, (0.0, 0.0, 0.0), 1.0), (29, (0.5, 0.5, 0.5), 1.0),
             (13, (0.0, 0.5, 0.25), 1.0), (13, (0.5, 0.0, 0.25), 1.0), (13, (0.0, 0.5, 0.75), 1.0), (13, (0.5, 0.0, 0.75), 1.0)]
    check("theta' conventional cell has 2 Cu + 4 Al", sorted(z for z, _, _ in sites) == [13, 13, 13, 13, 29, 29])
    return sites


def phases():
    th = theta_prime_sites()
    beta = beta_double_prime_sites()
    return [
        {"name": "Al", "cell": (4.0495,) * 3 + (90.0, 90.0, 90.0), "sites": [(13, p, 1.0) for p in FCC], "axes": CUBIC_AXES,
         "source": "apps preset Crystal.aluminum (mac4DSTEM 7d6fdbd1): fcc 4.0495, Z 13"},
        {"name": "Au", "cell": (4.0782,) * 3 + (90.0, 90.0, 90.0), "sites": [(79, p, 1.0) for p in FCC], "axes": CUBIC_AXES,
         "source": "apps preset Crystal.gold: fcc 4.0782, Z 79"},
        {"name": "Cu", "cell": (3.6149,) * 3 + (90.0, 90.0, 90.0), "sites": [(29, p, 1.0) for p in FCC], "axes": CUBIC_AXES,
         "source": "apps preset Crystal.copper: fcc 3.6149, Z 29"},
        {"name": "theta_prime_Al2Cu", "cell": (4.04, 4.04, 5.80, 90.0, 90.0, 90.0), "sites": th,
         "axes": [[0, 0, 1], [1, 0, 0], [0, 1, 1], [1, 1, 1], [1, 1, 2], [1, 2, 3]],
         "source": "no app or kit definition exists: body-centred tetragonal a 4.04, c 5.80, Cu 2a, Al 4d (fluorite-derived arrangement, sites written out); see theta_prime_sites"},
        {"name": "beta_double_prime_Mg5Si6", "cell": (15.16, 4.05, 6.74, 90.0, 105.3, 90.0), "sites": beta,
         "axes": CUBIC_AXES + [[0, 1, 0]],
         "source": "apps preset Crystal.betaDoublePrime (Andersen et al. 1998, Table 3 set 3, C2/m), sites in the apps' order; equals the CIF slice 2 pinned (cif_check)"},
        {"name": "Mg", "cell": (3.2094, 3.2094, 5.2108, 90.0, 90.0, 120.0), "sites": [(12, p, 1.0) for p in HCP],
         "axes": [[0, 0, 1], [1, 0, 0], [1, 1, 0], [2, 1, 0], [1, 1, 1]],
         "source": "apps preset Crystal.magnesium: hcp 3.2094 / 5.2108, Z 12; the kit's magnesium_P63mmc.cif cell"},
        {"name": "triclinic_test", "cell": (4.0, 5.0, 6.0, 80.0, 85.0, 95.0),
         "sites": [(13, (0.0, 0.0, 0.0), 1.0), (29, (0.31, 0.22, 0.13), 1.0)],
         "axes": [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 1, 1], [1, 1, 1]],
         "source": "the triclinic cell of registration-slice-1 (a 4, b 5, c 6, alpha 80, beta 85, gamma 95); the kit has NO triclinic CIF fixture, so two sites (Al at the origin, Cu at 0.31, 0.22, 0.13) are invented here"},
    ]


# ----------------------------------------------------------------------------------------------------------------
# The reference (V2): the definitions, in float64
# ----------------------------------------------------------------------------------------------------------------
def kit_frame(n):
    """The frame: e1 = normalise(n x r), e2 = n x e1, r = x-hat if |n_x| < 0.9 else y-hat."""
    r = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < BASIS_SWITCH else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(n, r)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return e1, e2


def cartesian_axis(A, uvw):
    v = uvw[0] * A[0] + uvw[1] * A[1] + uvw[2] * A[2]
    return v / np.linalg.norm(v)


def excited_entries(H, g, F2, n, e1, e2, kmax=K_MAX, slab=SLAB, width=WIDTH):
    """The pre-cap, pre-floor set of reference vectors. Returns arrays (all in the input order)."""
    s = g @ n
    q1, q2 = g @ e1, g @ e2
    qlen = np.sqrt(q1 * q1 + q2 * q2)
    w = max(width, 1e-9)
    I = F2 * np.exp(-((s / w) ** 2))
    sel = (np.abs(s) <= slab) & (qlen > Q_FLOOR) & (qlen <= kmax) & (I > 0)
    return sel, s, q1, q2, qlen, I


def library_entries(idx, H, q1, q2, qlen, I, floor=FLOOR, cap=CAP):
    """Floor, cap (strongest kept; ties: shorter |q| first), relative intensity, sort by |q| (ties: descending hkl)."""
    Imax = I[idx].max()
    keep = [i for i in idx if I[i] >= floor * Imax]
    n_floor = len(keep)
    tie_cap = False
    if len(keep) > cap:
        keep.sort(key=lambda i: (-I[i], qlen[i]))
        edge = I[keep[cap - 1]], I[keep[cap]]
        tie_cap = abs(edge[0] - edge[1]) <= 1e-9 * edge[0]
        keep = keep[:cap]
    keep.sort(key=lambda i: (qlen[i], tuple(-int(x) for x in H[i])))
    return keep, n_floor, tie_cap, Imax


def rot2(phi):
    c, s = math.cos(phi), math.sin(phi)
    return np.array([[c, -s], [s, c]])


def f(x):
    return float(x)


def fl(a):
    return [float(x) for x in a]


# ----------------------------------------------------------------------------------------------------------------
# Euler (SciPy) and the frame matrix B R_z(a)
# ----------------------------------------------------------------------------------------------------------------
def wrap_pi(a):
    """(-pi, pi]: SciPy's range is [-pi, pi]; -pi maps to +pi."""
    return math.pi if a <= -math.pi else float(a)


def euler_case(label, M):
    R = M.T                                         # crystal -> lab
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        ang = Rotation.from_matrix(R).as_euler("zxz")
    gimbal = any("imbal" in str(w.message) for w in wl)
    ang_w = [wrap_pi(a) for a in ang]
    Rb = Rotation.from_euler("zxz", ang_w).as_matrix()
    return {"label": label, "M": M.ravel().tolist(), "angles_raw": fl(ang), "angles": ang_w,
            "scipy_gimbal_warning": bool(gimbal), "crystal_to_lab": Rb.ravel().tolist(),
            "roundtrip_error": float(np.max(np.abs(Rb - R)))}


def euler_block(rng, orientation_matrices):
    cases = []
    rots = Rotation.random(1000, rng=rng).as_matrix()
    for i in range(1000):
        cases.append(euler_case(f"random-{i}", rots[i]))
    # gimbal and near-gimbal: Phi = 0, pi and 1e-9, 1e-10 away from each; 5e-9 (inside |sin Phi| <= 1e-8), then 2e-8, 1e-7,
    # 1e-6, 1e-5, 1e-3 (outside; a threshold of 1e-4 would treat 1e-5 and 1e-6 as locked, which the definition does not)
    pairs = [(0.7, -1.1), (-2.5, 2.0), (3.0, 3.0)]
    phis = [("0", 0.0), ("pi", math.pi), ("1e-9", 1e-9), ("1e-10", 1e-10), ("pi-1e-9", math.pi - 1e-9), ("pi-1e-10", math.pi - 1e-10),
            ("5e-9", 5e-9), ("2e-8", 2e-8), ("5e-8", 5e-8), ("9e-8", 9e-8), ("1e-7", 1e-7), ("1e-6", 1e-6), ("1e-5", 1e-5), ("1e-3", 1e-3),
            ("pi-5e-9", math.pi - 5e-9), ("pi-1e-5", math.pi - 1e-5)]
    gim = []
    for name, Phi in phis:
        for (p1, p2) in pairs:
            R = Rotation.from_euler("zxz", [p1, Phi, p2]).as_matrix()
            c = euler_case(f"gimbal-Phi={name}-phi1={p1}-phi2={p2}", R.T)
            c["input_angles"] = [p1, Phi, p2]
            c["Phi_label"] = name
            c["abs_sin_Phi"] = abs(math.sin(Phi))
            near_pi = abs(math.sin(Phi))                     # |sin Phi| distance to the nearest of 0, pi
            c["registered_rule_locks"] = bool(near_pi <= 1e-8)
            c["rules_agree"] = bool(c["registered_rule_locks"] == c["scipy_gimbal_warning"])
            gim.append(c)
    # py4DSTEM orientation matrices (parse_orientation, columns x', y', z'): [001] is Phi = 0 exactly
    realistic = [euler_case(lbl, np.array(M)) for lbl, M in orientation_matrices]
    return {"convention": "R = M^T (crystal -> lab), extrinsic z-x-z: Rotation.from_matrix(R).as_euler('zxz'), wrapped to (-pi, pi]; "
                          "crystal_to_lab = Rotation.from_euler('zxz', angles).as_matrix(); M is row-major (9 numbers)",
            "random": cases, "gimbal": gim, "py4dstem_orientation_matrices": realistic}


def rz_cases(rng, frames):
    """matrix(detectorBasis: B, inPlaneAngle: a) = B . R_z(a), R_z columns (c, -s, 0), (s, c, 0), (0, 0, 1)."""
    out = []
    bases = [(f"frame-{lbl}", B) for lbl, B in frames[:10]]
    rots = Rotation.random(10, rng=rng).as_matrix()
    bases += [(f"random-{i}", rots[i]) for i in range(10)]
    angles = rng.uniform(-math.pi, math.pi, size=20)
    angles[0], angles[1], angles[2] = 0.0, math.pi / 2, -math.pi
    for (lbl, B), a in zip(bases, angles):
        c, s = math.cos(a), math.sin(a)
        Rz = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])    # columns (c,-s,0), (s,c,0), (0,0,1)
        check("R_z columns are (c,-s,0),(s,c,0),(0,0,1)", np.allclose(Rz[:, 0], [c, -s, 0]) and np.allclose(Rz[:, 1], [s, c, 0]))
        out.append({"label": lbl, "B": np.asarray(B).ravel().tolist(), "angle": float(a), "matrix": (np.asarray(B) @ Rz).ravel().tolist()})
    return out


# ----------------------------------------------------------------------------------------------------------------
# Closed forms (NumPy from the written definitions)
# ----------------------------------------------------------------------------------------------------------------
def chance_match_fraction(V, rho, R):
    """p = min(1, V rho^2 / R^2), the matching's chance fraction. Defined here for V >= 0, rho >= 0, R > 0."""
    return min(1.0, V * rho * rho / (R * R))


def dist(a, b):
    dx, dy = a[0] - b[0], a[1] - b[1]
    return math.sqrt(dx * dx + dy * dy)


def nearest(u, refs, rho):
    """The closest reference within the radius (<=); ties: the lowest index. None when nothing is within rho."""
    best_i, best_d = None, None
    for i, r in enumerate(refs):
        d = dist(u, r)
        if best_d is None or d < best_d:
            best_i, best_d = i, d
    if best_d is not None and best_d <= rho:
        return best_i, best_d
    return None


def score_nearest_assignment(E, refs, rho):
    """Variant A: each experimental vector is assigned its nearest reference within rho; a reference's best distance is the
    least over the vectors assigned to it; the score is the mean over the unique references hit."""
    best = {}
    for e in E:
        r = nearest(e, refs, rho)
        if r is not None:
            i, d = r
            best[i] = min(best.get(i, d), d)
    return (math.fsum(best.values()) / len(best) if best else None), len(best), sorted(best)


def score_any_pair(E, refs, rho):
    """Variant B: a reference is hit when ANY experimental vector is within rho; its best distance is the least over them."""
    best = {}
    for i, r in enumerate(refs):
        ds = [dist(e, r) for e in E]
        d = min(ds) if ds else None
        if d is not None and d <= rho:
            best[i] = d
    return (math.fsum(best.values()) / len(best) if best else None), len(best), sorted(best)


def friedel_min(vs):
    m = None
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            d = dist((vs[i][0] + vs[j][0], vs[i][1] + vs[j][1]), (0.0, 0.0))
            m = d if m is None or d < m else m
    return m


def in_plane_steps(step_deg):
    """N = max(1, round(360/step)) (round half away from zero), theta_k = k 2 pi / N; one step theta = 0 for a step that is
    not finite, <= 0 or > 360."""
    if not math.isfinite(step_deg) or step_deg <= 0 or step_deg > 360:
        return [0.0]
    N = max(1, int(math.floor(360.0 / step_deg + 0.5)))
    return [k * 2 * math.pi / N for k in range(N)]


def low_index_axes():
    """u, v, w in -2..2 (u outer, then v, then w), not all zero, divided by their gcd, deduplicated up to sign, first
    representative kept."""
    seen, out = set(), []
    for u in range(-2, 3):
        for v in range(-2, 3):
            for w in range(-2, 3):
                if u == v == w == 0:
                    continue
                g = math.gcd(math.gcd(abs(u), abs(v)), abs(w))
                t = (u // g, v // g, w // g)
                neg = (-t[0], -t[1], -t[2])
                if t in seen or neg in seen:
                    continue
                seen.add(t)
                out.append([t[0], t[1], t[2]])
    return out


def closed_forms(rng):
    out = {}
    # chanceMatchFraction: formula cases (V, rho, R) incl. the clamp
    cm = []
    for V, rho, R in [(10, 0.02, 0.7), (48, 0.02, 0.7), (1, 0.05, 1.0), (100, 0.1, 0.5), (4000, 0.05, 1.0), (0, 0.02, 0.7),
                      (7, 0.0, 0.7), (48, 0.02, 1.6), (3, 0.015, 0.3), (1000, 0.02, 0.25), (25, 0.5, 0.5), (1, 1.0, 1.0),
                      (5, 0.3, 0.3)]:
        cm.append({"V": V, "rho": rho, "R": R, "p": chance_match_fraction(V, rho, R), "clamped": V * rho * rho / (R * R) > 1.0})
    for _ in range(12):
        V = int(rng.integers(1, 400)); rho = float(rng.uniform(0.005, 0.2)); R = float(rng.uniform(0.2, 1.6))
        cm.append({"V": V, "rho": rho, "R": R, "p": chance_match_fraction(V, rho, R), "clamped": V * rho * rho / (R * R) > 1.0})
    out["chanceMatchFraction"] = {
        "definition": "p = min(1, V rho^2 / R^2); cases have V >= 0, rho >= 0, R > 0",
        "cases": cm,
        "guards_UNREGISTERED": "the registration does not say what R <= 0, V <= 0 or rho <= 0 return; the cases above only use "
                               "V >= 0, rho >= 0, R > 0 (the formula gives 0 for V = 0 or rho = 0). The guard values of the "
                               "source are NOT pinned: the registration needs an amendment (proposed: R <= 0 or a non-finite "
                               "input returns 1, the safe 'cannot exclude chance')."}
    # nearest: 50 seeded cases
    nr = []
    for i in range(50):
        n_ref = int(rng.integers(5, 41))
        refs = rng.uniform(-0.8, 0.8, size=(n_ref, 2))
        rho = float([0.02, 0.05, 0.1][i % 3])
        k = int(rng.integers(0, n_ref))
        u = refs[k] + rng.normal(0, rho * 0.6, size=2) if i % 5 else rng.uniform(-0.8, 0.8, size=2)
        res = nearest(u, refs.tolist(), rho)
        ds = sorted(dist(u, r) for r in refs.tolist())
        tie = len(ds) > 1 and ds[1] - ds[0] <= 1e-12 * max(ds[0], 1e-300)
        edge = any(abs(d - rho) <= 1e-12 for d in ds)
        nr.append({"refs": refs.tolist(), "u": u.tolist(), "rho": rho, "index": None if res is None else res[0],
                   "distance": None if res is None else res[1], "tie": bool(tie), "edge": bool(edge),
                   "within_count": sum(1 for d in ds if d <= rho)})
    # the exact boundary (all numbers exactly representable): d = rho exactly -> a hit under <=
    nr.append({"refs": [[0.5, 0.0], [-0.5, 0.0]], "u": [0.75, 0.0], "rho": 0.25, "index": 0, "distance": 0.25, "tie": False,
               "edge": True, "within_count": 1, "label": "exact boundary: d == rho, a hit under <= (mutant 6 turns it into a miss)"})
    out["nearest"] = {"definition": "closest reference within rho (<=); ties: the lowest index (flag tie: skip them in a test that "
                                    "does not register the rule); edge: a reference at |d - rho| <= 1e-12 (skip unless the case is "
                                    "the exact boundary)", "cases": nr}
    check("no tie in the random nearest cases", not any(c["tie"] for c in nr))
    # score: 20 seeded cases (+ the exact boundary)
    sc = []
    for i in range(20):
        n_ref = int(rng.integers(6, 31))
        refs = rng.uniform(-0.8, 0.8, size=(n_ref, 2))
        rho = float([0.02, 0.03, 0.05][i % 3])
        E = []
        for k in rng.choice(n_ref, size=int(rng.integers(2, min(8, n_ref))), replace=False):
            E.append((refs[k] + rng.normal(0, rho * 0.5, size=2)).tolist())
            if i % 4 == 0:                           # a second vector on the same reference (best distance wins)
                E.append((refs[k] + rng.normal(0, rho * 0.8, size=2)).tolist())
        for _ in range(int(rng.integers(0, 4))):
            E.append(rng.uniform(-0.8, 0.8, size=2).tolist())          # an outlier
        if i % 7 == 6:
            E = [rng.uniform(-0.8, 0.8, size=2).tolist() for _ in range(3)]   # probably no hit
        a = score_nearest_assignment(E, refs.tolist(), rho)
        b = score_any_pair(E, refs.tolist(), rho)
        edge = any(abs(dist(e, r) - rho) <= 1e-12 for e in E for r in refs.tolist())
        sc.append({"refs": refs.tolist(), "E": E, "rho": rho, "score": a[0], "unique_hit": a[1], "hit_indices": a[2],
                   "variant_any_pair": {"score": b[0], "unique_hit": b[1]}, "variants_differ": bool(a[1] != b[1] or (a[0] != b[0] and not (a[0] is not None and b[0] is not None and abs(a[0] - b[0]) <= 1e-15 * abs(a[0])))),
                   "edge": bool(edge)})
    sc.append({"refs": [[0.5, 0.0], [-0.5, 0.0]], "E": [[0.75, 0.0], [-0.5, 0.0]], "rho": 0.25,
               "score": 0.125, "unique_hit": 2, "hit_indices": [0, 1], "variant_any_pair": {"score": 0.125, "unique_hit": 2},
               "variants_differ": False, "edge": True, "label": "exact boundary: d = 0.25 and 0.0, mean 0.125"})
    out["score"] = {"definition": "VARIANT A (pinned as 'score'): each experimental vector is assigned its nearest reference within "
                                  "rho; a reference's best distance is the least over its assigned vectors; the score is the mean over "
                                  "the unique references hit (null when none). The registration says 'the mean, over unique references "
                                  "hit, of their best distance': variant B (a reference is hit when any vector is within rho) is "
                                  "recorded as variant_any_pair; the cases where they differ are flagged variants_differ and a test "
                                  "must skip them unless the source's reading is confirmed.", "cases": sc}
    # containsFriedelPair: 20 cases
    fp = []
    for i in range(20):
        rho = float([0.02, 0.04][i % 2])
        n = int(rng.integers(3, 9))
        vs = rng.uniform(-0.8, 0.8, size=(n, 2)).tolist()
        kind = i % 4
        if kind in (0, 1):                           # a planted pair, noise below rho
            v = np.array(vs[0]); vs[1] = (-v + rng.normal(0, rho * 0.3, size=2)).tolist()
        elif kind == 2:                              # a pair that misses by ~1.5 rho
            v = np.array(vs[0]); ang = rng.uniform(0, 2 * math.pi); vs[1] = (-v + 1.5 * rho * np.array([math.cos(ang), math.sin(ang)])).tolist()
        m = friedel_min(vs)
        fp.append({"vectors": vs, "radius": rho, "contains": bool(m <= rho), "min_sum_norm": m, "edge": bool(abs(m - rho) <= 1e-12)})
    fp.append({"vectors": [[0.5, 0.0], [-0.25, 0.0], [0.1, 0.2]], "radius": 0.25, "contains": True, "min_sum_norm": 0.25, "edge": True,
               "label": "exact boundary: |v1 + v2| == radius, a pair under <="})
    out["containsFriedelPair"] = {
        "definition": "UNREGISTERED in section 4, so this is the reading proposed here: some pair i < j with |v_i + v_j| <= radius "
                      "(<=, as nearest). The pins record the minimum |v_i + v_j|. Flag edge: skip.", "cases": fp}
    # inPlaneSteps
    steps = [("2", 2.0), ("5", 5.0), ("7", 7.0), ("360", 360.0), ("0", 0.0), ("-1", -1.0), ("nan", float("nan")), ("400", 400.0),
             ("0.3", 0.3), ("11", 11.0), ("13", 13.0), ("inf", float("inf")), ("2.5", 2.5)]
    st = []
    for lbl, s in steps:
        a = in_plane_steps(s)
        if s == s and math.isfinite(s) and 0 < s <= 360:
            raw = 360.0 / s
            rnd_differs_trunc = int(math.floor(raw + 0.5)) != int(math.floor(raw))
        else:
            raw, rnd_differs_trunc = None, False
        st.append({"step": lbl, "count": len(a), "angles": a, "round_differs_from_truncation": bool(rnd_differs_trunc)})
    out["inPlaneSteps"] = {"definition": "N = max(1, round(360/step)), round half away from zero; theta_k = k 2 pi / N; [0] for a step "
                                         "that is not finite, <= 0 or > 360. Steps 11 and 13 are added so that rounding and truncation "
                                         "differ (360/11 = 32.7, 360/13 = 27.7; mutant 5).", "cases": st}
    # the 62 low-index axes
    ax = low_index_axes()
    ax_unreduced = []
    seen_u = set()
    for u in range(-2, 3):
        for v in range(-2, 3):
            for w in range(-2, 3):
                if (u, v, w) == (0, 0, 0) or (u, v, w) in seen_u or (-u, -v, -w) in seen_u:
                    continue
                seen_u.add((u, v, w))
                ax_unreduced.append([u, v, w])
    out["lowIndexZoneAxes"] = {"definition": "u, v, w in -2..2 (u outer, then v, then w), not all zero, / gcd, deduplicated up to "
                                             "sign, first representative kept", "count": len(ax), "axes": ax,
                               "REGISTRATION_INCONSISTENCY": "the registration says 62; the stated definition gives 49 (98 primitive "
                                                             "vectors in the box, / 2 for sign). 62 = 124 / 2 is the count of the "
                                                             "nonzero vectors deduplicated up to sign WITHOUT the gcd reduction "
                                                             "(then [200] and [100] are two 'axes'): listed as alternative_unreduced. "
                                                             "The apps-parity dump decides which the source computes.",
                               "alternative_unreduced_count": len(ax_unreduced), "alternative_unreduced": ax_unreduced}
    predict("the low-index zone axes number 62 (registration 3, 'Not settings')", 62, len(ax))
    predict("the unreduced sign-deduplicated vectors number 62 (the alternative reading)", 62, len(ax_unreduced))
    out["PhaseVectorResolution"] = {"status": "NOT PINNED: registration section 4 does not define PhaseVectorResolution, its pixel numbers or "
                                              "scaledToDetector, and the oracle is written from the registration, not from the Swift. "
                                              "Amendment needed: state the definitions (inputs, outputs, formulas) in section 4.",
                                    "cases": []}
    return out


# ----------------------------------------------------------------------------------------------------------------
# One phase
# ----------------------------------------------------------------------------------------------------------------
def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


def phase_block(ph, table, py, slice2):
    name, cell, sites, axes = ph["name"], ph["cell"], ph["sites"], ph["axes"]
    A, B, vol = lattice(cell)
    H, g = enumerate_box(A, B, K_MAX)
    order = order_reflections(H, g)
    H, g = H[order], g[order]
    gc = H.astype(float) @ B                                    # Cartesian g, Angstrom^-1, no 2 pi
    F = structure_factors(table, vol, H, g, sites)
    aF = np.abs(F)
    F2 = F.real ** 2 + F.imag ** 2
    kept = aF > CUT
    Hk, gk, gck, F2k, aFk = H[kept], g[kept], gc[kept], F2[kept], aF[kept]
    cut_edge = [list(map(int, Hk_)) for Hk_, a in zip(H, aF) if abs(a - CUT) <= 1e-7 * CUT]
    kmax_edge = [list(map(int, H_)) for H_, gg in zip(H, g) if abs(gg - K_MAX) <= 1e-12]
    print(f"[{name}] {len(H)} reflections within {K_MAX}, {int(kept.sum())} kept by |F| > {CUT}", flush=True)

    # ---- py4DSTEM's lattice and structure-factor table against this one ----
    p = py["phases"][name]
    lat_res = float(np.max(np.abs(np.array(p["lat_real"]) - A)))
    latinv_res = float(np.max(np.abs(np.array(p["lat_inv"]) - B)))
    check(f"{name}: py4DSTEM lat_real equals this lattice to 1e-12", lat_res <= 1e-12 * max(cell[:3]), str(lat_res))
    check(f"{name}: py4DSTEM lat_inv equals the reciprocal rows to 1e-12", latinv_res <= 1e-12, str(latinv_res))
    p_hkl = {tuple(h): i for i, h in enumerate(p["hkl"])}
    kit_set = {tuple(int(x) for x in h) for h in Hk}
    p_set = set(p_hkl)
    # py4DSTEM list: its own float32 cut; compare to the float64 kit cut
    missing = sorted(kit_set - p_set)                          # kept here, absent in py4DSTEM
    extra = sorted(p_set - kit_set - {(0, 0, 0)})               # in py4DSTEM, not kept here
    p_F2 = np.array(p["F2"])
    common = [tuple(int(x) for x in h) for h in Hk if tuple(int(x) for x in h) in p_hkl]
    f2_rel = [rel(p_F2[p_hkl[h]], F2k[i]) for i, h in enumerate(tuple(int(x) for x in hh) for hh in Hk) if h in p_hkl]
    p_g = np.array(p["g"])
    gdiff = max((float(np.max(np.abs(p_g[p_hkl[h]] - gck[i]))) for i, h in enumerate(tuple(int(x) for x in hh) for hh in Hk) if h in p_hkl), default=0.0)
    sfcmp = {"py4dstem_count": len(p_set), "kit_count": len(kit_set), "common": len(common),
             "kept_here_missing_in_py4dstem": [list(h) for h in missing], "in_py4dstem_not_kept_here": [list(h) for h in extra],
             "max_relative_F2_difference": max(f2_rel, default=0.0), "max_abs_g_difference": gdiff,
             "py4dstem_F_dtype": p["F_dtype"], "py4dstem_F2_dtype": p["F2_dtype"],
             "note": "py4DSTEM accumulates the structure factor in complex64, so its |F|^2 carries a float32 floor"}
    check(f"{name}: py4DSTEM's |F|^2 equals the float64 one to 1e-4 relative (the float32 floor)", sfcmp["max_relative_F2_difference"] <= 1e-4)
    # slice 2
    s2 = slice2.get(name)
    if s2 is not None:
        items = {(it[0], it[1], it[2]): it for it in s2["reflections"]["items"] if it[6] > CUT and it[3] <= 1.2}
        mine = {tuple(int(x) for x in h): i for i, h in enumerate(Hk) if gk[i] <= 1.2}
        cmp_f64 = [rel(F2k[mine[h]], it[6] ** 2) for h, it in items.items() if h in mine]
        cmp_p = [rel(p_F2[p_hkl[h]], it[6] ** 2) for h, it in items.items() if h in p_hkl]
        sfcmp["vs_slice2_pins_at_kMax_1_2"] = {
            "slice2_kept": len(items), "this_kept": len(mine), "hkl_only_in_slice2": sorted(map(list, set(items) - set(mine)))[:10],
            "hkl_only_here": sorted(map(list, set(mine) - set(items)))[:10],
            "max_relative_F2_difference_float64": max(cmp_f64, default=None), "max_relative_F2_difference_py4dstem_float32": max(cmp_p, default=None)}
        check(f"{name}: float64 |F|^2 equals slice 2's to 1e-11 relative (its sites are wrapped and CIF-rounded)", max(cmp_f64) <= 1e-11, str(max(cmp_f64)))
        check(f"{name}: the kept hkl sets at 1.2 equal slice 2's", set(items) == set(mine))

    zones = []
    obs_matrices = []
    frames = []
    worst = {"lambda0_vs_flat_relI": 0.0, "lambda0_vs_ref_relI": 0.0, "lambda0_vs_ref_q": 0.0, "frame_rotation": 0.0, "lambda0_vs_flat_px": 0.0}
    for uvw in axes:
        key = ",".join(str(x) for x in uvw)
        pa = p["axes"][key]
        O = np.array(pa["orientation"])
        n = cartesian_axis(A, uvw)
        e1, e2 = kit_frame(n)
        check(f"{name} [{key}]: (e1, e2, n) is right-handed orthonormal", np.allclose(np.cross(e1, e2), n, atol=1e-14) and abs(e1 @ e2) < 1e-15)
        check(f"{name} [{key}]: py4DSTEM's z column equals n to 1e-15", np.allclose(O[:, 2], n, atol=1e-15), str(O[:, 2] - n))
        xp, yp, zp = O[:, 0], O[:, 1], O[:, 2]
        M = np.array([[xp @ e1, xp @ e2], [yp @ e1, yp @ e2]])
        det = float(np.linalg.det(M))
        phi_py = math.atan2(M[1, 0], M[0, 0])
        hand = "rotation" if det > 0 else "reflection"
        check(f"{name} [{key}]: the 2x2 map is orthogonal", np.allclose(M @ M.T, np.eye(2), atol=1e-14))
        frames.append((f"{name}[{key}]", np.column_stack([e1, e2, n])))
        obs_matrices.append((f"py4dstem-{name}[{key}]", O))

        sel, s, q1, q2, qlen, I = excited_entries(Hk, gck, F2k, n, e1, e2)
        idx = [int(i) for i in np.where(sel)[0]]
        idx.sort(key=lambda i: (qlen[i], tuple(-int(x) for x in Hk[i])))
        slab_edge = [list(map(int, Hk[i])) for i in range(len(Hk)) if abs(abs(s[i]) - SLAB) <= 1e-12 and qlen[i] > Q_FLOOR]
        Imax = float(I[idx].max())

        # ---- V0 (py4DSTEM at lambda -> 0) and V1 (NumPy flat on py4DSTEM's own table) ----
        def py_set(rec):
            d = {}
            for j in range(len(rec["h"])):
                ql = math.hypot(rec["qx"][j], rec["qy"][j])
                if ql <= Q_FLOOR:
                    continue
                d[(rec["h"][j], rec["k"][j], rec["l"][j])] = (rec["qx"][j], rec["qy"][j], rec["qz"][j], rec["intensity"][j], ql)
            mx = max((v[3] for v in d.values()), default=1.0)
            return {h: (v[0], v[1], v[2], v[3], v[4], v[3] / mx) for h, v in d.items()}
        V0 = py_set(pa["lambda0"])
        # V1
        gp = np.array(p["g"]); Hp = np.array(p["hkl"]); F2p = np.array(p["F2"])
        s1 = gp @ zp
        px, py_ = gp @ xp, gp @ yp
        ql1 = np.sqrt(px * px + py_ * py_)
        I1 = F2p * np.exp(-((s1 / WIDTH) ** 2))
        sel1 = (np.abs(s1) <= SLAB) & (ql1 > Q_FLOOR) & (ql1 <= K_MAX) & (I1 > 0)
        mx1 = I1[sel1].max()
        V1 = {tuple(int(x) for x in Hp[i]): (px[i], py_[i], I1[i] / mx1) for i in np.where(sel1)[0]}
        check(f"{name} [{key}]: the lambda -> 0 py4DSTEM set equals the flat set on py4DSTEM's own table", set(V0) == set(V1),
              f"{sorted(set(V0) ^ set(V1))}")
        r_I = max((rel(V0[h][5], V1[h][2]) for h in V0), default=0.0)
        r_p = max((max(abs(V0[h][0] - V1[h][0]), abs(V0[h][1] - V1[h][1])) for h in V0), default=0.0)
        worst["lambda0_vs_flat_relI"] = max(worst["lambda0_vs_flat_relI"], r_I)
        worst["lambda0_vs_flat_px"] = max(worst["lambda0_vs_flat_px"], r_p)

        # ---- V2 against V0 ----
        ref_set = {tuple(int(x) for x in Hk[i]) for i in idx}
        only_ref = sorted(ref_set - set(V0))
        only_py = sorted(set(V0) - ref_set)
        both = ref_set & set(V0)
        d_I = max((rel(V0[h][5], I[ii] / Imax) for ii in idx for h in [tuple(int(x) for x in Hk[ii])] if h in V0), default=0.0)
        d_q = max((abs(V0[h][4] - qlen[ii]) for ii in idx for h in [tuple(int(x) for x in Hk[ii])] if h in V0), default=0.0)
        Rm = rot2(phi_py) if det > 0 else rot2(phi_py) @ np.diag([1.0, -1.0])
        d_rot = 0.0
        for ii in idx:
            h = tuple(int(x) for x in Hk[ii])
            if h in V0:
                pred = Rm @ np.array([q1[ii], q2[ii]])
                d_rot = max(d_rot, float(max(abs(pred[0] - V0[h][0]), abs(pred[1] - V0[h][1]))))
        worst["lambda0_vs_ref_relI"] = max(worst["lambda0_vs_ref_relI"], d_I)
        worst["lambda0_vs_ref_q"] = max(worst["lambda0_vs_ref_q"], d_q)
        worst["frame_rotation"] = max(worst["frame_rotation"], d_rot)

        # ---- the pinned entries ----
        entries = []
        for ii in idx:
            h = tuple(int(x) for x in Hk[ii])
            v = V0.get(h)
            entries.append([h[0], h[1], h[2], f(gck[ii][0]), f(gck[ii][1]), f(gck[ii][2]), f(s[ii]), f(q1[ii]), f(q2[ii]), f(qlen[ii]),
                            f(F2k[ii]), f(I[ii]), f(I[ii] / Imax),
                            None if v is None else v[0], None if v is None else v[1], None if v is None else v[4], None if v is None else v[5]])
        keep, n_floor, tie_cap, _ = library_entries(idx, Hk, q1, q2, qlen, I)
        library = [[int(Hk[i][0]), int(Hk[i][1]), int(Hk[i][2]), f(q1[i]), f(q2[i]), f(qlen[i]), f(I[i] / Imax)] for i in keep]
        floor_edge = [list(map(int, Hk[i])) for i in idx if abs(I[i] / Imax - FLOOR) <= 1e-9]

        # ---- V200 ----
        V2c = py_set(pa["kv200"])
        enter = sorted(set(V2c) - ref_set)
        leave = sorted(ref_set - set(V2c))
        d200 = max((abs(V2c[h][5] - I[ii] / Imax) for ii in idx for h in [tuple(int(x) for x in Hk[ii])] if h in V2c), default=0.0)
        zones.append({
            "uvw": list(uvw), "n": fl(n), "e1": fl(e1), "e2": fl(e2),
            "py4dstem": {"orientation_matrix_columns_xyz": O.tolist(), "x": fl(xp), "y": fl(yp), "z": fl(zp),
                         "wavelength_lambda0": pa["wavelength_lambda0"], "wavelength_200kV": pa["wavelength_200kV"]},
            "frame_map": {"definition": "(g.x', g.y') = M (g.e1, g.e2) with M = [[x'.e1, x'.e2], [y'.e1, y'.e2]]; for det +1, M = R(phi_py) = [[c,-s],[s,c]]; "
                                        "for det -1, M = R(phi_py) diag(1,-1)",
                          "M": M.ravel().tolist(), "phi_py_rad": phi_py, "phi_py_deg": math.degrees(phi_py), "det": det, "handedness": hand},
            "counts": {"pre_cap_pre_floor": len(idx), "after_floor": n_floor, "after_cap": len(keep), "py4dstem_lambda0": len(V0),
                       "only_reference": [list(h) for h in only_ref], "only_py4dstem": [list(h) for h in only_py],
                       "max_relative_I_vs_py4dstem": d_I, "max_abs_q_vs_py4dstem": d_q, "max_abs_rotated_q_vs_py4dstem": d_rot,
                       "lambda0_vs_flat_numpy_on_py4dstem_table": {"max_relative_I": r_I, "max_abs_px": r_p}},
            "flags": {"slab_edge_hkl": slab_edge, "floor_edge_hkl": floor_edge, "cap_tie": bool(tie_cap), "cut_edge_hkl": cut_edge,
                      "kmax_edge_hkl": kmax_edge},
            "entry_fields": ["h", "k", "l", "gx", "gy", "gz", "s", "q1", "q2", "qlen", "F2", "I_excited", "I_relative",
                             "py4dstem_x", "py4dstem_y", "py4dstem_qlen", "py4dstem_I_relative"],
            "entries": entries,
            "library_fields": ["h", "k", "l", "q1", "q2", "qlen", "I_relative"],
            "library": library,
            "kv200": {"count_curved": len(V2c), "count_flat": len(ref_set), "entering": [list(h) for h in enter],
                      "leaving": [list(h) for h in leave], "max_abs_dI_relative_common": d200,
                      "definition": "entering = in py4DSTEM's 200 kV curved set but not in the flat slab; leaving = the reverse"},
        })
    Hlist = [[int(a) for a in h] for h in Hk]
    return {"name": name, "source": ph["source"],
            "cell": dict(zip(["a", "b", "c", "alpha", "beta", "gamma"], [float(x) for x in cell])), "volume": vol,
            "lat_real": A.tolist(), "lat_inv": B.tolist(),
            "sites": [{"z": int(z), "frac": [float(x) for x in pos], "occupancy": float(occ)} for (z, pos, occ) in sites],
            "kept_reflections": len(Hlist), "reflections_within_kMax": int(len(H)),
            "structure_factor_check": sfcmp, "lattice_check": {"max_abs_lat_real_difference_vs_py4dstem": lat_res, "max_abs_lat_inv_difference_vs_py4dstem": latinv_res},
            "zone_axes": zones, "worst": worst}, frames, obs_matrices


# ----------------------------------------------------------------------------------------------------------------
def run_py4dstem(request):
    check("the py4DSTEM interpreter exists", os.path.isfile(PY4D_PYTHON), PY4D_PYTHON)
    r = subprocess.run([PY4D_PYTHON, "-I", "-B", PY4D_HELPER], input=json.dumps(request).encode(), capture_output=True, timeout=1200)
    if r.returncode != 0:
        sys.stderr.write(r.stderr.decode())
        raise AssertionError(f"SELF-CHECK FAILED: the py4DSTEM helper exited {r.returncode}")
    return json.loads(r.stdout)


def main():
    table = load_table()
    PH = phases()
    for ph in PH:
        ph["_sites_json"] = [{"z": z, "frac": list(p), "occ": o} for (z, p, o) in ph["sites"]]
    request = {"k_max": K_MAX, "slab": SLAB, "width": WIDTH,
               "phases": [{"name": ph["name"], "cell": list(ph["cell"]), "sites": ph["_sites_json"], "axes": ph["axes"]} for ph in PH]}
    py = run_py4dstem(request)
    s2 = {s["name"]: s for s in json.load(open(SLICE2_PINS, encoding="utf-8"))["structures"]}
    slice2 = {"Al": s2["Al"], "Au": s2["Au"], "Mg": s2["Mg"], "beta_double_prime_Mg5Si6": s2["beta_double_prime"]}

    blocks, frames, mats = [], [], []
    for ph in PH:
        b, fr, mt = phase_block(ph, table, py, slice2)
        blocks.append(b)
        frames += fr
        mats += mt
    beta = next(ph for ph in PH if ph["name"] == "beta_double_prime_Mg5Si6")
    cif_check = beta_cif_check(beta["sites"])
    next(b for b in blocks if b["name"] == "beta_double_prime_Mg5Si6")["cif_check"] = cif_check

    # ---- aggregate residuals ----
    agg = {k: max(b["worst"][k] for b in blocks) for k in blocks[0]["worst"]}
    # the registered G2 tolerances
    predict("py4DSTEM lambda -> 0 residual vs NumPy-flat on its own table <= 1e-12 relative (registration 4/5)", True,
            bool(agg["lambda0_vs_flat_relI"] <= 1e-12))
    predict("py4DSTEM relative intensity vs the float64 reference <= 1e-9 relative (G2 as registered)", True,
            bool(agg["lambda0_vs_ref_relI"] <= 1e-9))
    predict("the kit's q rotated by phi_py equals py4DSTEM's within 1e-12 (G2)", True, bool(agg["frame_rotation"] <= 1e-12))
    predict("|q| vs py4DSTEM within 1e-12 (G2)", True, bool(agg["lambda0_vs_ref_q"] <= 1e-12))
    for b in blocks:
        for z in b["zone_axes"]:
            c = z["counts"]
            predict(f"{b['name']} {z['uvw']}: the hkl set equals py4DSTEM's", 0, len(c["only_reference"]) + len(c["only_py4dstem"]))
        predict(f"{b['name']}: py4DSTEM's tiling misses no kept hkl", 0, len(b["structure_factor_check"]["kept_here_missing_in_py4dstem"]))
        predict(f"{b['name']}: py4DSTEM keeps no hkl the float64 cut drops", 0, len(b["structure_factor_check"]["in_py4dstem_not_kept_here"]))

    # ---- Euler and R_z ----
    rng = np.random.default_rng(SEED)
    eul = euler_block(rng, mats)
    rzs = rz_cases(rng, frames)
    check("1000 random Euler cases", len(eul["random"]) == 1000)
    worst_rt = max(c["roundtrip_error"] for c in eul["random"] + eul["py4dstem_orientation_matrices"])
    check("every random and py4DSTEM Euler round trip (angles -> matrix) returns R to 1e-12", worst_rt <= 1e-12, str(worst_rt))
    # a locked case (SciPy sets the third angle to 0) returns R only up to the discarded Phi: the round trip is not exact there
    eul["gimbal_roundtrip_worst_by_Phi"] = {}
    for c in eul["gimbal"]:
        d = eul["gimbal_roundtrip_worst_by_Phi"]
        d[c["Phi_label"]] = max(d.get(c["Phi_label"], 0.0), c["roundtrip_error"])
    scipy_gimbal = {c["Phi_label"]: [] for c in eul["gimbal"]}
    for c in eul["gimbal"]:
        scipy_gimbal[c["Phi_label"]].append(c["scipy_gimbal_warning"])
    eul["scipy_gimbal_warning_by_Phi"] = {k: v for k, v in scipy_gimbal.items()}
    disagree = sorted({c["Phi_label"] for c in eul["gimbal"] if not c["rules_agree"]})
    eul["gimbal_rule_note"] = ("SciPy 1.18.1 locks when the middle angle is below 1e-7 (a scan of 400 random (phi1, phi2) each: 100% lock at "
                               "9.9e-8, 50% at exactly 1e-7, 0% at 1.01e-7), not at |sin Phi| <= 1e-8 as registration 4 says. The cases whose "
                               "Phi_label is in rule_disagreement have rules_agree false: a Swift port that follows the registered 1e-8 gives "
                               "phi2 != 0 there while SciPy gives 0 (the same rotation, different angles). Skip them or amend the registration "
                               "to SciPy's rule; Phi = 1e-5 and 1e-6 are unlocked in both, so mutant 13 (threshold 1e-4) is still caught.")
    eul["rule_disagreement"] = disagree
    predict("the registered gimbal rule (|sin Phi| <= 1e-8) agrees with SciPy 1.18.1 on every pinned gimbal case", [], disagree)
    eul["roundtrip_worst"] = worst_rt

    cf = closed_forms(rng)

    # ---- the header and the record ----
    out = {
        "header": {
            "purpose": "slice 4 (vector matching and phase mapping, DiffractoOrientation 0.3.0): oracle pins, registration docs/registration-slice-4.md section 6",
            "generator": "tools/slice-4/pins.py + tools/slice-4/pins-py4dstem.py",
            "commands": ["References/python/orix/bin/python -I -B DiffractoKit/tools/slice-4/pins.py",
                         "(started by pins.py) References/python/py4dstem/bin/python -I -B DiffractoKit/tools/slice-4/pins-py4dstem.py < request.json"],
            "versions": {"pins.py": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                                     "pymatgen": dist_version("pymatgen"), "pymatgen_core": dist_version("pymatgen-core")},
                         "pins-py4dstem.py": py["versions"]},
            "py4dstem_commit": py["provenance"]["lock"], "py4dstem_tree_head": py["provenance"]["tree_head"],
            "py4dstem_provenance": py["provenance"], "py4dstem_call": py["call"],
            "scattering_table_sha256": TABLE_SHA256,
            "seed": SEED,
            "registered_defaults": {"kMaxInvAngstrom": K_MAX, "excitationSlabInvAngstrom": SLAB, "excitationWidthInvAngstrom": WIDTH,
                                    "minimumIntensityFraction": FLOOR, "maximumVectorsPerEntry": CAP, "inPlaneStepDeg": IN_PLANE_STEP,
                                    "F_cut": CUT, "q_floor": Q_FLOOR, "basis_switch": BASIS_SWITCH},
            "oracles": {"V2_reference": "NumPy float64 from the registration's definitions (these pins' 'entries' and 'library')",
                        "V0": "py4DSTEM Crystal.generate_diffraction_pattern at a huge accelerating voltage (lambda ~ 1e-26 A), sigma = w/sqrt(2), "
                              "tol_excitation_error_mult = slab/sigma, tol_intensity 0: the 'py4dstem_*' fields of 'entries'",
                        "V1": "NumPy flat excitation on py4DSTEM's own hkl, g and float32 |F|^2: the lambda -> 0 residual",
                        "V200": "py4DSTEM at 200 kV, curved Ewald: the 'kv200' records"},
            "files_read": ["Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json", "References/py4dstem-scattering/scattering_factors.txt",
                           "References/crystal-structures/beta_double_prime_Mg5Si6.cif (the beta'' site check only)"],
        },
        "phases": blocks,
        "frame": {"definition": "e1 = normalise(n x r), e2 = n x e1, r = x-hat if |n_x| < 0.9 else y-hat, n the Cartesian zone axis "
                                "from the real lattice rows; q = (g.e1, g.e2); py4DSTEM: parse_orientation(zone_axis_lattice), p = (g.x', g.y')"},
        "measured": {"worst_over_all_phases_and_axes": agg,
                     "float32_note": "py4DSTEM's structure factors are complex64: lambda0_vs_ref_relI is that floor, "
                                     "lambda0_vs_flat_relI is the exp/threshold arithmetic alone (both share the float32 |F|^2)"},
        "euler": eul,
        "frame_matrix_cases": {"definition": "matrix = B R_z(a), R_z = [[c, s, 0], [-s, c, 0], [0, 0, 1]] (columns (c,-s,0), (s,c,0), (0,0,1)), "
                                             "row-major 9 numbers; the first 10 B are the columns (e1, e2, n) of the frame of the first 10 "
                                             "phase/axis pairs (an assumption about the basis layout: the product does not depend on it), "
                                             "then 10 seeded random rotations",
                               "cases": rzs},
        "closed_forms": cf,
        "predictions": [{"name": n, "predicted": p_, "measured": m_, "equal": p_ == m_} for n, p_, m_ in PREDICTIONS],
        "self_checks": len(CHECKS),
    }

    os.makedirs(OUT_DIR, exist_ok=True)
    dump = lambda o: json.dumps(o, separators=(",", ":"), sort_keys=True, allow_nan=False)
    with open(OUT, "w", encoding="utf-8") as fh:
        keys = sorted(out)
        fh.write("{\n")
        for i, k in enumerate(keys):
            if k == "phases":
                fh.write('"phases":[\n' + ",\n".join(dump(b) for b in out["phases"]) + "\n]")
            elif k == "euler":
                e = out["euler"]
                fh.write('"euler":{' + ",\n".join(f'"{kk}":' + (("[\n" + ",\n".join(dump(c) for c in e[kk]) + "\n]") if isinstance(e[kk], list) else dump(e[kk]))
                                                  for kk in sorted(e)) + "}")
            else:
                fh.write(f'"{k}":' + dump(out[k]))
            fh.write(",\n" if i < len(keys) - 1 else "\n")
        fh.write("}\n")
    size = os.path.getsize(OUT)
    json.load(open(OUT, encoding="utf-8"))
    print(f"\nwrote {OUT}  ({size / 1e6:.2f} MB, {len(blocks)} phases, {sum(len(b['zone_axes']) for b in blocks)} zone axes)")
    print(f"{len(CHECKS)} self-checks ran, all passed (a failure would have stopped the run)")
    print("worst residuals:", json.dumps(agg))
    bad = [(n, p_, m_) for n, p_, m_ in PREDICTIONS if p_ != m_]
    print(f"{len(PREDICTIONS)} predictions, {len(PREDICTIONS) - len(bad)} equal, {len(bad)} differ")
    for n, p_, m_ in bad:
        print(f"!!! PREDICTION DIFFERS: {n}: registered {p_}, measured {m_}")
    if bad:
        sys.exit(3)


if __name__ == "__main__":
    main()
