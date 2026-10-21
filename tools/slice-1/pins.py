"""DiffractoKit slice 1: the pin generator (registration docs/registration-slice-1.md, sections 3, 4, 5.1, 5.3, amendment 8).

Run, from any directory, with the pinned oracle:

    ../References/python/orix/bin/python -I -B tools/slice-1/pins.py

It writes Tests/DiffractoCoreTests/Fixtures/slice-1-pins.json BEFORE any Swift exists. pymatgen 2026.9.24 (with
pymatgen-core 2026.9.23, spglib 2.8.0, NumPy 2.5.4, Python 3.14.7) is the oracle for the cell, the expanded sites, the
space-group number, the point-group symbol, the symmetry operations, the (h, k, l) set with |g| <= k_max and d.
The Laue orbit, the canonical form, the absence test and the families are NOT pymatgen functions: they are written
here, in NumPy, from the definitions of registration section 3 (3.5, 3.6, 3.8, 3.9), not from any Swift.
Agreement with pymatgen is agreement, not truth; truth is registration section 5 (closed forms are checked below too).

File:line pins into the pymatgen-core clone (References/pymatgen-core-ae9a834/src/pymatgen/):
  core/lattice.py       Lattice.from_parameters (vesta=True branch) 389-445, metric_tensor 206, volume 182,
                        reciprocal_lattice_crystallographic 539, d_hkl 258, get_points_in_sphere 1390
  io/cif.py             CifParser.parse_structures 1298, _get_structure 970 (site merge atol = site_tolerance 1e-4 at 1037)
  symmetry/analyzer.py  SpacegroupAnalyzer.get_symmetry_operations 338 (rotations and translations in fractional
                        coordinates; _get_symmetry rounds the translations to Fraction.limit_denominator(1000)),
                        get_space_group_number 166, get_point_group_symbol 194
  core/operations.py    SymmOp.as_xyz_str 462 (via util/string.py transformation_to_string 325)
"""
import datetime
import hashlib
import itertools
import json
import math
import os
import re
import sys
import warnings
from importlib.metadata import version as dist_version

import numpy as np
import spglib

import pymatgen.core
from pymatgen.core import Lattice, Structure
from pymatgen.core.operations import SymmOp
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# ----------------------------------------------------------------------------------------------------------------
# Paths (relative to this script, so it runs from any directory)
# ----------------------------------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURES = os.path.join(KIT, "Tests", "DiffractoCoreTests", "Fixtures")
CRYSTALS = os.path.join(KIT, "References", "crystal-structures")
OUT = os.path.join(FIXTURES, "slice-1-pins.json")

# ----------------------------------------------------------------------------------------------------------------
# The registered constants
# ----------------------------------------------------------------------------------------------------------------
K_MAX = 1.2                 # Angstrom^-1 (registration 4)
EXTRA_K_MAX = (1.6, 2.0)    # beta'' only, counts only
FAMILY_TOL = 1e-9           # registration 3.9
TIE_TOL = 1e-12             # registration 3.8: equal |g| within 1e-12 relative
SYMPREC = 1e-4              # registration 2 / 4
ORBIT_ITEMS = 24            # reflections per structure that carry their "orbit" in the pins (canonical is on every item)
ABSENCE_TOL = 1e-9          # registration 3.8: (h k l) . t an integer within 1e-9
HEX_AXES_CLASSES = ("-3", "-3m", "6/m", "6/mmm")   # registration 3.6: four Miller-Bravais indices

# The 32 -> 11 table. Keys are spglib/pymatgen point-group symbols.
LAUE_OF_POINT_GROUP = {
    "1": "-1", "-1": "-1",
    "2": "2/m", "m": "2/m", "2/m": "2/m",
    "222": "mmm", "mm2": "mmm", "mmm": "mmm",
    "4": "4/m", "-4": "4/m", "4/m": "4/m",
    "422": "4/mmm", "4mm": "4/mmm", "-42m": "4/mmm", "-4m2": "4/mmm", "4/mmm": "4/mmm",
    "3": "-3", "-3": "-3",
    "32": "-3m", "3m": "-3m", "-3m": "-3m",
    "6": "6/m", "-6": "6/m", "6/m": "6/m",
    "622": "6/mmm", "6mm": "6/mmm", "-6m2": "6/mmm", "-62m": "6/mmm", "6/mmm": "6/mmm",
    "23": "m-3", "m-3": "m-3",
    "432": "m-3m", "-43m": "m-3m", "m-3m": "m-3m",
}

# Cell contents the registration (5.3) states, per fixture (the cell content the file's formula states).
EXPECTED_SITES = {
    "beta_double_prime": 22, "beta_double_prime_P1": 22,
    "gold_P1": 4, "gold_centring_ops": 4, "silicon_Fd-3m": 8, "magnesium_P63mmc": 2,
    "rocksalt_Fm-3m_192ops": 8, "alpha_quartz_P3221": 9, "calcite_R-3c_hex": 30, "fluorapatite_P63m": 42,
    "pyrite_Pa-3": 12, "ws2_2H": 6, "cdi2_P-3m1": 3, "tetragonal_a_equals_c": 2, "ideal_hcp_P1": 2,
    # extras beyond the 5.3 list
    "magnesium_P63mmc_24ops": 2, "tetragonal_a_equals_c_4ops": 5, "alpha_quartz_P3221_as_seeded": 12, "pyrite_Pa-3_as_seeded": 8,
}
# Fixtures in the order of registration 5.3 (monoclinic first), then the extras.
FIXTURE_ORDER = [
    "gold_P1", "gold_centring_ops", "silicon_Fd-3m", "magnesium_P63mmc", "rocksalt_Fm-3m_192ops",
    "alpha_quartz_P3221", "calcite_R-3c_hex", "fluorapatite_P63m", "pyrite_Pa-3", "ws2_2H", "cdi2_P-3m1",
    "tetragonal_a_equals_c", "ideal_hcp_P1",
    "magnesium_P63mmc_24ops", "tetragonal_a_equals_c_4ops", "alpha_quartz_P3221_as_seeded", "pyrite_Pa-3_as_seeded",
]

# Registration 5.1: the seven parameter cells (monoclinic first).
PARAMETER_CELLS = [
    ("monoclinic", (15.16, 4.05, 6.74, 90.0, 105.3, 90.0)),
    ("cubic", (4.0, 4.0, 4.0, 90.0, 90.0, 90.0)),
    ("tetragonal", (3.0, 3.0, 5.0, 90.0, 90.0, 90.0)),
    ("orthorhombic", (3.0, 4.0, 5.0, 90.0, 90.0, 90.0)),
    ("hexagonal", (3.0, 3.0, 5.0, 90.0, 90.0, 120.0)),
    ("rhombohedral", (4.0, 4.0, 4.0, 70.0, 70.0, 70.0)),
    ("triclinic", (4.0, 5.0, 6.0, 80.0, 85.0, 95.0)),
]

# The files whose operation loop is complete: pymatgen's detected set must equal what the file lists.
FILE_OPS_MUST_EQUAL_DETECTED = {
    "beta_double_prime", "rocksalt_Fm-3m_192ops", "alpha_quartz_P3221", "calcite_R-3c_hex", "pyrite_Pa-3",
    "magnesium_P63mmc_24ops",
}

CHECKS = []   # (name, ok) for the summary; every check is also an assert


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    assert ok, f"SELF-CHECK FAILED: {name} {detail}"


# ----------------------------------------------------------------------------------------------------------------
# Lattice (registration 3.1)
# ----------------------------------------------------------------------------------------------------------------
def lattice_block(a, b, c, alpha, beta, gamma):
    """pymatgen's vesta=True lattice (a along x, b in xy) and everything the pin file records about it."""
    L = Lattice.from_parameters(a, b, c, alpha, beta, gamma, vesta=True)
    A = np.array(L.matrix)
    G = L.metric_tensor
    V = L.volume
    B = np.array(L.reciprocal_lattice_crystallographic.matrix)
    tag = f"cell({a},{b},{c},{alpha},{beta},{gamma})"

    # a along x, b in the xy plane
    check(f"{tag}: a along x, b in xy", abs(A[0, 1]) == 0 and abs(A[0, 2]) == 0 and abs(A[1, 2]) == 0)
    # reciprocal rows: a_i . b_j = delta_ij
    check(f"{tag}: a_i . b_j = delta_ij", np.allclose(A @ B.T, np.eye(3), rtol=0, atol=1e-12))
    # registration 3.1: B = G^-1 A and G = A A^T
    check(f"{tag}: B = G^-1 A", np.allclose(B, np.linalg.inv(G) @ A, rtol=0, atol=1e-12))
    # volume: |det A| equals the closed form to 1e-12 relative
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    v_closed = a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    check(f"{tag}: volume closed form", abs(V - v_closed) <= 1e-12 * v_closed, f"{V} {v_closed}")
    # registration 3.1's third row
    sg = math.sin(math.radians(gamma))
    a3 = np.array([c * cb, -c * (cb * cg - ca) / sg, V / (a * b * sg)])
    check(f"{tag}: third lattice row closed form", np.allclose(A[2], a3, rtol=0, atol=1e-12 * max(a, b, c)))
    return L, A, G, V, B


def g_vectors(H, B):
    """g(hkl) = h b1 + k b2 + l b3 for an (M, 3) integer array; returns the vectors and |g|."""
    gv = np.asarray(H, dtype=float) @ B
    return gv, np.linalg.norm(gv, axis=1)


def hkl_range_table(L, B, lo=-3, hi=3):
    """Every h, k, l in lo..hi except 0,0,0 with |g| and d (registration 5.1)."""
    H = np.array([h for h in itertools.product(range(lo, hi + 1), repeat=3) if h != (0, 0, 0)], dtype=np.int64)
    _, g = g_vectors(H, B)
    d = 1.0 / g
    for i in range(len(H)):
        ref = L.d_hkl(tuple(int(x) for x in H[i]))
        check("table d equals Lattice.d_hkl to 1e-12 relative", abs(d[i] - ref) <= 1e-12 * ref, f"{H[i]} {d[i]} {ref}")
    return [{"hkl": [int(x) for x in H[i]], "g": float(g[i]), "d": float(d[i])} for i in range(len(H))]


# ----------------------------------------------------------------------------------------------------------------
# Truth for 5.1: closed forms (independent of pymatgen)
# ----------------------------------------------------------------------------------------------------------------
def truth_inv_d2(system, params, h, k, l):
    a, b, c, al, be, ga = params
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    sa, sb, sg = (math.sin(math.radians(x)) for x in (al, be, ga))
    if system == "cubic":
        return (h * h + k * k + l * l) / (a * a)
    if system == "tetragonal":
        return (h * h + k * k) / (a * a) + l * l / (c * c)
    if system == "orthorhombic":
        return h * h / (a * a) + k * k / (b * b) + l * l / (c * c)
    if system == "hexagonal":
        return (4.0 / 3.0) * (h * h + h * k + k * k) / (a * a) + l * l / (c * c)
    if system == "rhombohedral":
        num = (h * h + k * k + l * l) * sa * sa + 2 * (h * k + k * l + h * l) * (ca * ca - ca)
        den = a * a * (1 - 3 * ca * ca + 2 * ca ** 3)
        return num / den
    if system == "monoclinic":
        return (h * h / (a * a) + l * l / (c * c) - 2 * h * l * cb / (a * c)) / (sb * sb) + k * k / (b * b)
    if system == "triclinic":
        return truth_triclinic(params, h, k, l)
    raise ValueError(system)


def truth_triclinic(params, h, k, l):
    a, b, c, al, be, ga = params
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    sa, sb, sg = (math.sin(math.radians(x)) for x in (al, be, ga))
    V = a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    s11 = b * b * c * c * sa * sa
    s22 = a * a * c * c * sb * sb
    s33 = a * a * b * b * sg * sg
    s12 = a * b * c * c * (ca * cb - cg)
    s23 = a * a * b * c * (cb * cg - ca)
    s13 = a * b * b * c * (ca * cg - cb)
    return (s11 * h * h + s22 * k * k + s33 * l * l + 2 * s12 * h * k + 2 * s23 * k * l + 2 * s13 * h * l) / (V * V)


def check_truth_5_1(system, params, table):
    worst = 0.0
    worst_tri = 0.0
    for e in table:
        h, k, l = e["hkl"]
        d_formula = 1.0 / math.sqrt(truth_inv_d2(system, params, h, k, l))
        d_tri = 1.0 / math.sqrt(truth_triclinic(params, h, k, l))
        worst = max(worst, abs(e["d"] - d_formula) / d_formula)
        worst_tri = max(worst_tri, abs(e["d"] - d_tri) / d_tri)
    check(f"5.1 {system}: closed form equals the metric-tensor d to 1e-12 relative", worst <= 1e-12, f"{worst}")
    check(f"5.1 {system}: the general triclinic formula equals it too", worst_tri <= 1e-12, f"{worst_tri}")
    return worst, worst_tri


def check_monoclinic_anchors(table):
    # macTEM S2c hand anchors, registration 5.1 (5-6 significant digits as written there)
    want = {(2, 0, 0): 7.31135, (0, 0, 1): 6.50112, (0, 2, 0): 2.025, (1, 1, 0): 3.90306,
            (2, 0, 1): 4.32457, (2, 0, -1): 5.65554}
    by = {tuple(e["hkl"]): e["d"] for e in table}
    for hkl, d in want.items():
        check(f"5.1 monoclinic anchor d{hkl} = {d}", abs(by[hkl] - d) <= 6e-6, f"{by[hkl]}")


# ----------------------------------------------------------------------------------------------------------------
# Symmetry: operations, Laue group, Laue class (registration 3.3, 3.5)
# ----------------------------------------------------------------------------------------------------------------
I3 = np.eye(3, dtype=np.int64)


def matrix_order(R):
    M = I3.copy()
    for n in range(1, 25):
        M = M @ R
        if np.array_equal(M, I3):
            return n
    raise AssertionError("no finite order")


def snap_translation(t):
    """Mod 1 into [0, 1); snapped to the nearest k/24 when within 1e-9 (every crystallographic translation is a
    multiple of 1/24; pymatgen already rounds spglib's to Fraction.limit_denominator(1000))."""
    out = []
    for x in t:
        k24 = x * 24.0
        k = round(k24)
        check("translation is a multiple of 1/24 within 1e-9", abs(k24 - k) <= 24e-9, f"{x}")
        out.append((k % 24) / 24.0)
    return out


def operations_of(structure):
    """Rotations (int 3x3), snapped translations, xyz strings, as pymatgen's SpacegroupAnalyzer reports them."""
    sga = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    ops = sga.get_symmetry_operations(cartesian=False)
    rots, trans, xyz = [], [], []
    for op in ops:
        R = np.array(op.rotation_matrix)
        check("rotation part is an integer matrix (|R - round(R)| < 1e-8)", np.max(np.abs(R - np.round(R))) < 1e-8)
        R = np.round(R).astype(np.int64)
        t = snap_translation(np.array(op.translation_vector, dtype=float))
        rots.append(R)
        trans.append(t)
        xyz.append(SymmOp.from_rotation_and_translation(R, t).as_xyz_str().replace(" ", ""))
    return sga, np.array(rots), np.array(trans), xyz


def op_key(R, t):
    return (tuple(int(x) for x in R.ravel()), tuple(int(round(x * 24)) % 24 for x in t))


def check_ops_closed_and_distinct(name, Rs, ts):
    keys = [op_key(R, t) for R, t in zip(Rs, ts)]
    check(f"{name}: operations are distinct", len(set(keys)) == len(keys))
    kset = set(keys)
    ok = True
    for R1, t1 in zip(Rs, ts):
        for R2, t2 in zip(Rs, ts):
            # (R1, t1) after (R2, t2): p -> R1 (R2 p + t2) + t1
            if op_key(R1 @ R2, R1 @ np.array(t2) + np.array(t1)) not in kset:
                ok = False
                break
        if not ok:
            break
    check(f"{name}: operations are closed under composition", ok)


def check_ops_permute_sites(name, Rs, ts, Z, occ, frac):
    """Every operation maps the site set onto itself (Z and occupancy preserved), positions mod 1 within 1e-4."""
    n = len(Z)
    for R, t in zip(Rs, ts):
        P = (frac @ R.T + t)
        used = set()
        for i in range(n):
            d = P[i] - frac
            d -= np.round(d)
            hit = np.where((np.max(np.abs(d), axis=1) < 1e-4) & (Z == Z[i]) & (np.abs(occ - occ[i]) < 1e-9))[0]
            assert len(hit) >= 1, f"{name}: operation {op_key(R, t)} does not map site {i} onto a site"
            used.add(int(hit[0]))
        assert len(used) == n, f"{name}: operation {op_key(R, t)} is not a permutation of the sites"
    check(f"{name}: every operation maps the site set onto itself", True)


def laue_rotations(Rs):
    s = {tuple(int(x) for x in R.ravel()) for R in Rs} | {tuple(int(-x) for x in R.ravel()) for R in Rs}
    return np.array(sorted(s), dtype=np.int64).reshape(-1, 3, 3)


def check_group_closed(name, Rs):
    s = {tuple(int(x) for x in R.ravel()) for R in Rs}
    ok = all(tuple(int(x) for x in (A @ B).ravel()) in s for A in Rs for B in Rs)
    check(f"{name}: Laue rotation set is a group (closed)", ok)


def laue_class_from_group(Rs, proper_only=True):
    """Registration 3.5: by the order and the axis orders present. 'A 4-fold / a 6-fold' is read as a PROPER rotation
    (det +1) of that order; read over every element (improper ones included) the table fails, because -3 has order 6
    (see the report): proper_only=False reproduces the literal reading for the diagnostic."""
    n = len(Rs)
    dets = [int(round(np.linalg.det(R))) for R in Rs]
    orders = [matrix_order(R) for R in Rs]
    sel = [o for o, d in zip(orders, dets) if (d == 1 or not proper_only)]
    has4 = 4 in sel
    has6 = 6 in sel
    if n == 2:
        return "-1"
    if n == 4:
        return "2/m"
    if n == 8:
        return "4/m" if has4 else "mmm"
    if n == 16:
        return "4/mmm"
    if n == 6:
        return "-3"
    if n == 12:
        return "6/m" if has6 else "-3m"
    if n == 24:
        return "6/mmm" if has6 else "m-3"
    if n == 48:
        return "m-3m"
    return f"?order{n}"


# ----------------------------------------------------------------------------------------------------------------
# Reflections (registration 3.6, 3.8, 3.9)
# ----------------------------------------------------------------------------------------------------------------
def enumerate_hkl(L, A, B, kmax):
    """pymatgen's get_points_in_sphere on the crystallographic reciprocal lattice, cross-checked against a brute-force
    box |h| <= ceil(kmax |a_i|) (the registered per-axis bound)."""
    rec = L.reciprocal_lattice_crystallographic
    pts = rec.get_points_in_sphere([[0.0, 0.0, 0.0]], [0.0, 0.0, 0.0], kmax)
    py = set()
    py_dist = {}
    for fc, dist, idx, img in pts:
        h = tuple(int(round(x)) for x in fc)
        assert max(abs(x - round(x)) for x in fc) < 1e-9
        if h != (0, 0, 0):
            py.add(h)
            py_dist[h] = float(dist)
    lens = np.linalg.norm(A, axis=1)
    nmax = [int(math.ceil(kmax * x)) for x in lens]
    grid = np.array(list(itertools.product(*[range(-n, n + 1) for n in nmax])), dtype=np.int64)
    grid = grid[np.any(grid != 0, axis=1)]
    _, g = g_vectors(grid, B)
    bf = {tuple(int(x) for x in row) for row in grid[g <= kmax]}
    check(f"hkl set |g| <= {kmax}: get_points_in_sphere equals the brute-force box", py == bf,
          f"only pymatgen: {sorted(py - bf)[:5]} only box: {sorted(bf - py)[:5]}")
    worst_g = max(abs(py_dist[h] - float(np.linalg.norm(np.array(h, dtype=float) @ B))) / py_dist[h] for h in py)
    check(f"|g| at k_max {kmax}: pymatgen's get_points_in_sphere distance equals |h b1 + k b2 + l b3| to 1e-12 relative",
          worst_g <= 1e-12, f"{worst_g}")
    return sorted(py)


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
        cluster = sorted(idx[i:j].tolist(), key=lambda m: tuple(int(x) for x in H[m]), reverse=True)
        out.extend(cluster)
        i = j
    return out


def canonical_key(m):
    return (sum(1 for x in m if x >= 0), m)


class OrbitCache:
    def __init__(self, laue_R):
        self.R = laue_R
        self.cache = {}

    def orbit(self, hkl):
        o = self.cache.get(hkl)
        if o is None:
            imgs = np.einsum("i,nij->nj", np.array(hkl, dtype=np.int64), self.R)
            o = sorted({tuple(int(x) for x in row) for row in imgs})
            self.cache[hkl] = o
        return o

    def canonical(self, hkl):
        return max(self.orbit(hkl), key=canonical_key)


def allowed_flags(H, Rs, ts):
    """Registration 3.8: absent iff some operation has (hkl) R == (hkl) and (hkl) . t not an integer within 1e-9."""
    HR = np.einsum("mi,nij->mnj", H, Rs)
    fixed = np.all(HR == H[:, None, :], axis=2)
    dot = H.astype(float) @ np.asarray(ts, dtype=float).T
    nonint = np.abs(dot - np.round(dot)) > ABSENCE_TOL
    return ~np.any(fixed & nonint, axis=1)


def orbit_name(c, hex_axes):
    h, k, l = c
    if hex_axes:
        return "{%d %d %d %d}" % (h, k, -(h + k), l)
    return "{%d %d %d}" % (h, k, l)


def build_families(seq, cache, hex_axes):
    """seq: [(hkl tuple, g)] already in sweep order. Registration 3.9."""
    fams = []
    cur = []
    g0 = None

    def close():
        members = sorted((m for m, _ in cur), reverse=True)
        orbits = sorted({cache.canonical(m) for m in members}, reverse=True)
        g = min(gg for _, gg in cur)
        fams.append({
            "g": float(g), "d": float(1.0 / g),
            "members": [list(m) for m in members],
            "multiplicity": len(members),
            "orbits": [list(o) for o in orbits],
            "name": " + ".join(orbit_name(o, hex_axes) for o in orbits),
        })

    for hkl, g in seq:
        if cur and g > g0 * (1 + FAMILY_TOL):
            close()
            cur = []
        if not cur:
            g0 = g
        cur.append((hkl, g))
    if cur:
        close()
    return fams


def reflections_for(L, A, B, kmax, Rs, ts, laue_R, hex_axes, cache, detail):
    hkl = enumerate_hkl(L, A, B, kmax)
    H = np.array(hkl, dtype=np.int64)
    _, g = g_vectors(H, B)
    d = 1.0 / g
    worst_d = max(abs(d[i] - L.d_hkl(tuple(int(x) for x in H[i]))) / d[i] for i in range(len(H)))
    check(f"every reflection at k_max {kmax}: d = 1/|g| equals Lattice.d_hkl to 1e-12 relative", worst_d <= 1e-12, f"{worst_d}")
    order = order_reflections(H, g)
    Ho = H[order]
    go = g[order]
    allowed = allowed_flags(Ho, Rs, ts)

    # near-coincidence diagnostic: relative gaps between successive distinct |g| clusters in (1e-12, 1e-6]
    rel_gaps = []
    prev = go[0]
    for x in go[1:]:
        if x > prev * (1 + TIE_TOL):
            rel_gaps.append((x - prev) / prev)
            prev = x
    near = sorted(r for r in rel_gaps if r <= 1e-6)

    seq_all = [(tuple(int(x) for x in Ho[i]), float(go[i])) for i in range(len(Ho))]
    seq_allowed = [s for s, a in zip(seq_all, allowed) if a]
    fams = build_families(seq_allowed, cache, hex_axes)
    fams_all = build_families(seq_all, cache, hex_axes)
    res = {"count": len(Ho), "allowed_count": int(allowed.sum()), "families_n": len(fams), "families": fams,
           "families_all_n": len(fams_all), "families_all": fams_all, "near_coincident": near}
    if detail:
        items = []
        for i, (m, gg) in enumerate(seq_all):
            item = {"hkl": list(m), "g": gg, "d": float(1.0 / gg), "allowed": bool(allowed[i])}
            if i < ORBIT_ITEMS:       # the orbit is kept for the first ORBIT_ITEMS reflections only (file size)
                item["orbit"] = [list(o) for o in cache.orbit(m)]
            item["canonical"] = list(cache.canonical(m))
            items.append(item)
        res["items"] = items
    return res


# ----------------------------------------------------------------------------------------------------------------
# Truth for 5.3 (group theory and cell contents, independent of pymatgen's numbers)
# ----------------------------------------------------------------------------------------------------------------
def sorted_abs(h):
    return tuple(sorted((abs(x) for x in h), reverse=True))


def check_truth_5_3(name, entry):
    items = entry["reflections"]["items"]
    allowed = {tuple(i["hkl"]): i["allowed"] for i in items}

    def iff(label, rule):
        bad = [h for h in allowed if allowed[h] != rule(*h)]
        check(f"5.3 {name}: {label}", not bad, f"{bad[:5]}")

    if name in ("beta_double_prime", "beta_double_prime_P1"):
        iff("allowed iff h + k even (C-centring)", lambda h, k, l: (h + k) % 2 == 0)
        check("5.3 beta'': Laue 2/m, 8 operations, 22 sites, space group 12",
              entry["laue_class"] == "2/m" and len(entry["operations"]) == 8 and entry["space_group_number"] == 12)
        check("5.3 beta'': 1448 allowed reflections and 419 families at 1.2 (registered prediction)",
              entry["reflections"]["allowed_count"] == 1448 and len(entry["families"]) == 419,
              f"{entry['reflections']['allowed_count']} {len(entry['families'])}")
        want = [("{2 0 0}", 2), ("{0 0 1}", 2), ("{2 0 -1}", 2), ("{2 0 1}", 2), ("{1 1 0}", 4), ("{4 0 0}", 2),
                ("{4 0 -1}", 2), ("{1 1 -1}", 4), ("{2 0 -2}", 2), ("{0 0 2}", 2)]
        got = [(f["name"], f["multiplicity"]) for f in entry["families"][:10]]
        check("5.3 beta'': the first ten families (registered)", got == want, f"{got}")
        if name == "beta_double_prime":
            ex = entry["extra_k_max"]
            check("5.3 beta'': extra k_max entries present", set(ex) == {"1.6", "2.0"})
    if name in ("gold_P1", "gold_centring_ops", "rocksalt_Fm-3m_192ops"):
        iff("allowed iff h, k, l all of one parity (F-centring)", lambda h, k, l: h % 2 == k % 2 == l % 2)
        check(f"5.3 {name}: space group 225, Laue m-3m, 192 operations",
              entry["space_group_number"] == 225 and entry["laue_class"] == "m-3m" and len(entry["operations"]) == 192)
    if name in ("gold_P1", "gold_centring_ops"):
        names = [f["name"] for f in entry["families"]]
        check("5.3 gold: 8 families {1 1 1} {2 0 0} {2 2 0} {3 1 1} {2 2 2} {4 0 0} {3 3 1} {4 2 0} (registered)",
              names == ["{1 1 1}", "{2 0 0}", "{2 2 0}", "{3 1 1}", "{2 2 2}", "{4 0 0}", "{3 3 1}", "{4 2 0}"], f"{names}")
    if name == "silicon_Fd-3m":
        absent = {(2, 0, 0), (4, 2, 0), (6, 0, 0)}
        present = {(2, 2, 2), (4, 4, 2)}
        for h, a in allowed.items():
            k = sorted_abs(h)
            if k in absent:
                assert not a, (name, h)
            if k in present:
                assert a, (name, h)
        check("5.3 silicon: {2 0 0} {4 2 0} {6 0 0} absent, {2 2 2} {4 4 2} allowed by the operations", True)
        check("5.3 silicon: 8 sites, m-3m, 227, 12 families (registered)",
              len(entry["sites"]) == 8 and entry["laue_class"] == "m-3m" and entry["space_group_number"] == 227
              and len(entry["families"]) == 12)
    if name in ("magnesium_P63mmc", "magnesium_P63mmc_24ops", "ws2_2H", "ideal_hcp_P1", "fluorapatite_P63m"):
        ok = all(allowed[h] == (h[2] % 2 == 0) for h in allowed if h[0] == 0 and h[1] == 0)
        check(f"5.3 {name}: (0 0 l) allowed iff l even (6_3 screw)", ok)
    if name in ("magnesium_P63mmc", "magnesium_P63mmc_24ops"):
        check("5.3 magnesium: 6/mmm, 194, {1 0 -1 0} 6, {0 0 0 2} 2, {1 0 -1 1} 12, {1 1 -2 0} 6",
              entry["laue_class"] == "6/mmm" and entry["space_group_number"] == 194
              and {f["name"]: f["multiplicity"] for f in entry["families"]}.items()
              >= {"{1 0 -1 0}": 6, "{0 0 0 2}": 2, "{1 0 -1 1}": 12, "{1 1 -2 0}": 6}.items())
    if name == "alpha_quartz_P3221":
        ok = all(allowed[h] == (h[2] % 3 == 0) for h in allowed if h[0] == 0 and h[1] == 0)
        check("5.3 quartz: (0 0 l) allowed iff l = 3n (3_2 screw)", ok)
        check("5.3 quartz: 9 sites, Laue -3m (not 6/mmm), 154", entry["laue_class"] == "-3m" and entry["space_group_number"] == 154)
    if name == "calcite_R-3c_hex":
        ok = all(not a for h, a in allowed.items() if (-h[0] + h[1] + h[2]) % 3 != 0)
        check("5.3 calcite: -h + k + l not a multiple of 3 is absent (R-centring)", ok)
        check("5.3 calcite: 30 sites, -3m, 167, 36 operations",
              entry["laue_class"] == "-3m" and entry["space_group_number"] == 167 and len(entry["operations"]) == 36)
    if name == "fluorapatite_P63m":
        check("5.3 fluorapatite: 42 sites, 6/m (not 6/mmm), 176", entry["laue_class"] == "6/m" and entry["space_group_number"] == 176)
        orb = {tuple(f["orbits"][0]): f for f in entry["families"]}
        # {2 1 0}-type: the two orbits {2 1 -3 0} and {1 2 -3 0} are distinct under 6/m (24 would be one orbit)
        fam = [f for f in entry["families"] if f["name"] == "{2 1 -3 0} + {1 2 -3 0}"]
        check("5.3 fluorapatite: {2 1 0} and {1 2 0} are two orbits (6/m), one ring of 12", len(fam) == 1 and fam[0]["multiplicity"] == 12)
    if name == "pyrite_Pa-3":
        iff("Pa-3: absent iff (h = 0, k odd) or (k = 0, l odd) or (l = 0, h odd)",
            lambda h, k, l: not ((h == 0 and k % 2) or (k == 0 and l % 2) or (l == 0 and h % 2)))
        check("5.3 pyrite: 12 sites, m-3 (not m-3m), 205", entry["laue_class"] == "m-3" and entry["space_group_number"] == 205)
        fa = entry["families_all"]
        # the two orbits of the (2 1 0) ring under m-3 have canonical forms (2 1 0) and (2 0 1) by 3.6's key (the cyclic
        # permutation carries (1 2 0) to (2 0 1)); the registration's prose calls them {2 1 0} and {1 2 0}
        fam_all = [n for n in fa["first_names"] if n == "{2 1 0} + {2 0 1}"]
        check("5.3 pyrite: the (2 1 0) ring is two orbits of 12 on the full list", len(fam_all) == 1, f"{fa['first_names']}")
    if name == "ws2_2H":
        check("5.3 WS2: 6 sites, 6/mmm, 194", entry["laue_class"] == "6/mmm" and entry["space_group_number"] == 194)
    if name == "cdi2_P-3m1":
        check("5.3 CdI2: 3 sites, -3m, 164, no absences", entry["laue_class"] == "-3m" and entry["space_group_number"] == 164
              and entry["reflections"]["allowed_count"] == entry["reflections"]["count"])
    if name == "ideal_hcp_P1":
        fam = [f for f in entry["families"] if f["orbits"] == [[3, 0, 2], [0, 0, 6]]]
        check("5.3 ideal hcp: {3 0 -3 2} + {0 0 0 6} is one ring of 14", len(fam) == 1 and fam[0]["multiplicity"] == 14
              and fam[0]["name"] == "{3 0 -3 2} + {0 0 0 6}")
        check("5.3 ideal hcp: 6/mmm, 194", entry["laue_class"] == "6/mmm" and entry["space_group_number"] == 194)
    if name == "tetragonal_a_equals_c":
        check("5.3 tetragonal a = c: the Laue class is read from the atoms (mmm), not refused, not m-3m",
              entry["laue_class"] not in ("m-3m", "m-3"))


# ----------------------------------------------------------------------------------------------------------------
# One structure
# ----------------------------------------------------------------------------------------------------------------
def sha256_of(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def declared_number(text):
    m = re.search(r"^\s*_(?:symmetry_Int_Tables_number|space_group_IT_number)\s+(\d+)\s*$", text, re.M)
    return int(m.group(1)) if m else None


def raw_cell_from_cif(parser):
    blocks = parser.as_dict()
    assert len(blocks) == 1, "one data block expected"
    data = next(iter(blocks.values()))

    def f(tag):
        s = str(data[tag])
        return float(re.sub(r"\(.*\)", "", s))

    return tuple(f(t) for t in ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                                "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))


def wrap_unit(x):
    v = x % 1.0
    if v >= 1.0 - 1e-12:
        v = 0.0
    return float(v)


def sites_of(structure):
    Z, occ, frac = [], [], []
    for site in structure:
        assert site.is_ordered, "disordered site"
        (el, amt), = site.species.items()
        Z.append(int(el.Z))
        occ.append(float(amt))
        frac.append([wrap_unit(x) for x in site.frac_coords])
    return np.array(Z, dtype=np.int64), np.array(occ), np.array(frac)


def max_site_shift(ZA, fA, ZB, fB):
    """Largest periodic distance (fractional, per component) from each site of A to its partner of B (same Z, unique)."""
    assert len(ZA) == len(ZB)
    used = set()
    worst = 0.0
    for i in range(len(ZA)):
        d = fB - fA[i]
        d -= np.round(d)
        dist = np.max(np.abs(d), axis=1)
        cand = [j for j in np.argsort(dist) if ZB[j] == ZA[i] and j not in used]
        j = cand[0]
        used.add(int(j))
        worst = max(worst, float(dist[j]))
    return worst


def structure_entry(name, cell, Z, occ, frac, sym_structure, cif_file, cif_sha, cif_text, declared, derived_from,
                    pmg_warnings, extra_kmax, written=None, file_ops=None):
    a, b, c, al, be, ga = cell
    L, A, G, V, B = lattice_block(a, b, c, al, be, ga)
    expected = EXPECTED_SITES.get(name)
    if expected is not None:
        check(f"{name}: {expected} sites (the cell content the registration states)", len(Z) == expected, f"{len(Z)}")

    sga, Rs, ts, xyz = operations_of(sym_structure)
    sg_number = int(sga.get_space_group_number())
    pg = sga.get_point_group_symbol()
    check(f"{name}: point-group symbol {pg!r} is in the 32 -> 11 table", pg in LAUE_OF_POINT_GROUP)
    laue_table = LAUE_OF_POINT_GROUP[pg]
    if declared is not None:
        check(f"{name}: declared space group {declared} equals spglib's {sg_number}", declared == sg_number)
    check_ops_closed_and_distinct(name, Rs, ts)
    check_ops_permute_sites(name, Rs, ts, Z, occ, frac)

    file_ops_info = None
    if file_ops is not None:
        fkeys = set()
        for op in file_ops:
            R = np.round(np.array(op.rotation_matrix)).astype(np.int64)
            fkeys.add(op_key(R, np.array(op.translation_vector, dtype=float)))
        dkeys = {op_key(R, t) for R, t in zip(Rs, ts)}
        file_ops_info = {"count": len(fkeys), "equals_detected": fkeys == dkeys, "subset_of_detected": fkeys <= dkeys}
        if name in FILE_OPS_MUST_EQUAL_DETECTED:
            check(f"{name}: the file's own {len(fkeys)} operations equal the {len(dkeys)} spglib detects (hand-completed loop is right)",
                  fkeys == dkeys, f"{len(fkeys)} vs {len(dkeys)}")

    LR = laue_rotations(Rs)
    check_group_closed(name, LR)
    laue_group = laue_class_from_group(LR)
    check(f"{name}: Laue class from the rotation set (3.5) equals the table's ({laue_group} vs {laue_table})",
          laue_group == laue_table)
    literal_reading = laue_class_from_group(LR, proper_only=False)   # diagnostic only
    laue = laue_table
    hex_axes = laue in HEX_AXES_CLASSES

    cache = OrbitCache(LR)
    main = reflections_for(L, A, B, K_MAX, Rs, ts, LR, hex_axes, cache, detail=True)

    # self-checks on orbits: closed (orbit of any member is the orbit), canonical constant on an orbit, distinct between
    canon_to_orbit = {}
    for it in main["items"]:
        hkl = tuple(it["hkl"])
        orb = cache.orbit(hkl)
        can = cache.canonical(hkl)
        for m in orb:
            assert cache.orbit(m) == orb, f"{name}: orbit of {m} differs from the orbit of {hkl}"
            assert cache.canonical(m) == can
        prev = canon_to_orbit.setdefault(can, orb)
        assert prev == orb, f"{name}: canonical form {can} names two different orbits"
    check(f"{name}: orbit is closed, canonical form is constant on an orbit and differs between orbits", True)

    entry = {
        "name": name,
        "cif_file": cif_file,
        "cif_sha256": cif_sha,
    }
    if cif_text is not None:
        entry["cif_text"] = cif_text
    if derived_from is not None:
        entry["derived_from"] = derived_from
    if written is not None:
        entry["_written"] = written
    entry.update({
        "cell": {"a": a, "b": b, "c": c, "alpha": al, "beta": be, "gamma": ga},
        "matrix": A.tolist(), "metric": G.tolist(), "volume": float(V), "reciprocal": B.tolist(),
        "sites": [{"z": int(Z[i]), "frac": [float(x) for x in frac[i]], "occupancy": float(occ[i])}
                  for i in range(len(Z))],
        "declared_space_group_number": declared,
        "space_group_number": sg_number,
        "point_group": pg,
        "laue_class": laue,
        "operations": [{"rotation": Rs[i].tolist(), "translation": [float(x) for x in ts[i]], "xyz": xyz[i]} for i in range(len(Rs))],
        "laue_rotations": LR.tolist(),
        "reflections": {"k_max": K_MAX, "count": main["count"], "allowed_count": main["allowed_count"],
                        "items": main["items"]},
        "families": main["families"],
        "families_all": {"count": main["families_all_n"], "first_names": [f["name"] for f in main["families_all"][:12]]},
        "pymatgen_warnings": pmg_warnings,
    })
    if file_ops_info is not None:
        entry["file_operations"] = file_ops_info
    extras = {}
    for km in extra_kmax:
        r = reflections_for(L, A, B, km, Rs, ts, LR, hex_axes, cache, detail=False)
        extras[repr(km)] = {"count": r["count"], "allowed_count": r["allowed_count"], "families": r["families_n"]}
    if extras:
        entry["extra_k_max"] = extras

    wr = entry.pop("_written", None)
    if wr is not None:
        Zw, occw, fracw, st_w = wr
        shift = max_site_shift(Z, frac, Zw, fracw)
        entry["site_snap_max"] = shift
        if shift > 1e-10:
            # pymatgen rounds a coordinate within 1e-4 of 1/3 or 2/3 to the exact fraction before expanding
            # (io/cif.py _frac_tolerance, around line 610); the sites as written are pinned beside it
            entry["sites_as_written"] = [{"z": int(Zw[i]), "frac": [float(x) for x in fracw[i]], "occupancy": float(occw[i])}
                                         for i in range(len(Zw))]
            sga_w = SpacegroupAnalyzer(st_w, symprec=SYMPREC)
            check(f"{name}: the sites as written (no 1/3, 2/3 snap) give the same space group and operation count",
                  int(sga_w.get_space_group_number()) == sg_number and len(sga_w.get_symmetry_operations()) == len(Rs))
    check_truth_5_3(name, entry)

    # stdout
    print(f"\n=== {name}")
    print(f"space group {sg_number}  point group {pg}  Laue class {laue}  sites {len(Z)}  operations {len(Rs)}"
          f"  (declared {declared}; Laue rotations {len(LR)}; 3.5 read over all elements would say {literal_reading})")
    print(f"reflections at {K_MAX}: {main['count']} / allowed {main['allowed_count']}   families {main['families_n']}"
          f"   (families on the full list: {main['families_all_n']})")
    for f in main["families"][:10]:
        print(f"   {f['name']:<24s} d = {f['d']:.5f}  mult {f['multiplicity']}")
    for km, r in extras.items():
        print(f"at {km}: reflections {r['count']} / allowed {r['allowed_count']}   families {r['families']}")
    if main["near_coincident"]:
        print(f"   NOTE successive distinct |g| closer than 1e-6 relative: {['%.3e' % x for x in main['near_coincident'][:6]]}")
    if pmg_warnings:
        print(f"   pymatgen warnings: {len(pmg_warnings)}")
    return entry, {"name": name, "sg": sg_number, "pg": pg, "laue": laue, "sites": len(Z), "ops": len(Rs),
                   "count": main["count"], "allowed": main["allowed_count"], "families": main["families_n"],
                   "literal_laue": literal_reading, "declared": declared}


def parse_cif(path):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        parser = CifParser(path)
        structs = parser.parse_structures(primitive=False)
    msgs = sorted({str(x.message).strip().splitlines()[0] for x in w})
    assert len(structs) == 1, f"{path}: {len(structs)} structures"
    return parser, structs[0], msgs


def fixture_entry(name, path, cif_file, include_text, extra_kmax=()):
    text = open(path, encoding="utf-8").read()
    parser, st, msgs = parse_cif(path)
    cell = raw_cell_from_cif(parser)
    # pymatgen's own lattice parameters agree with the CIF's numbers to 1e-12 relative
    for x, y in zip(cell, st.lattice.parameters):
        check(f"{name}: cell parameter {x} equals pymatgen's {y} to 1e-12 relative", abs(x - y) <= 1e-12 * max(abs(x), 1.0))
    Z, occ, frac = sites_of(st)
    # the same file with pymatgen's 1/3, 2/3 snapping switched off (frac_tolerance = 0): the sites as written
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st_w = CifParser(path, frac_tolerance=0.0).parse_structures(primitive=False)[0]
    Zw, occw, fracw = sites_of(st_w)
    file_ops = parser.symmetry_operations if (re.search(r"^\s*_(symmetry_equiv_pos_as_xyz|space_group_symop_operation_xyz)", text, re.M)) else None
    if include_text:
        check(f"{name}: cif_text reproduces the file byte for byte (sha256)", hashlib.sha256(text.encode("utf-8")).hexdigest() == sha256_of(path))
    return structure_entry(name, cell, Z, occ, frac, st, cif_file, sha256_of(path), text if include_text else None,
                           declared_number(text), None, msgs, extra_kmax, written=(Zw, occw, fracw, st_w), file_ops=file_ops)


def beta_p1_entry(path):
    """The beta'' cell as a list of the 22 expanded sites given to pymatgen as P1: the symmetry is detected, not read."""
    parser, st, msgs = parse_cif(path)
    cell = raw_cell_from_cif(parser)
    Z, occ, frac = sites_of(st)
    L = Lattice.from_parameters(*cell, vesta=True)
    p1 = Structure(L, [int(z) for z in Z], frac, coords_are_cartesian=False)
    check("beta_double_prime_P1: the P1 structure carries 22 sites and no symmetry from a file", len(p1) == 22)
    return structure_entry("beta_double_prime_P1", cell, Z, occ, frac, p1, None, None, None, None,
                           "beta_double_prime", [], ())


# ----------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------
def main():
    versions = {
        "python": sys.version,
        "pymatgen": dist_version("pymatgen"),
        "pymatgen_core": dist_version("pymatgen-core"),
        "spglib": spglib.__version__,
        "numpy": np.__version__,
    }
    check("pymatgen.core.__version__ equals the pymatgen-core distribution", pymatgen.core.__version__ == versions["pymatgen_core"])
    print("oracle:", json.dumps(versions))

    out = {
        "versions": versions,
        "generated_at": datetime.datetime.now(datetime.UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "definitions": {
            "k_max": K_MAX, "family_tolerance": FAMILY_TOL, "symprec": SYMPREC,
            "laue_orbit": "hkl @ R over the rotation parts, plus Friedel",
            "canonical": "max by (count of non-negative indices, lexicographic (h,k,l))",
            "tie_tolerance": TIE_TOL, "absence_tolerance": ABSENCE_TOL,
            "translations": "mod 1 into [0,1), snapped to the nearest k/24 (asserted within 1e-9 of it)",
            "laue_class_rule": "registration 3.5 with 'a 4-fold / a 6-fold' read as a proper rotation (det +1) of that order",
            "family_name_order": "orbits by descending lexicographic canonical form, joined by ' + '",
            "reflection_order": "increasing |g|; ties (within 1e-12 relative of the cluster's first) by descending (h,k,l)",
        },
    }

    # ---- the seven parameter cells (lattice only) ----
    cells = []
    print("\n--- parameter cells (registration 5.1): closed forms against the metric tensor ---")
    for system, params in PARAMETER_CELLS:
        L, A, G, V, B = lattice_block(*params)
        table = hkl_range_table(L, B)
        w1, w2 = check_truth_5_1(system, params, table)
        if system == "monoclinic":
            check_monoclinic_anchors(table)
        a, b, c, al, be, ga = params
        cells.append({"name": system, "a": a, "b": b, "c": c, "alpha": al, "beta": be, "gamma": ga,
                      "matrix": A.tolist(), "metric": G.tolist(), "volume": float(V), "reciprocal": B.tolist(), "d": table})
        print(f"{system:<13s} V = {V:.9f}  max |d - formula|/d = {w1:.2e}  (triclinic formula {w2:.2e})")
    out["cells"] = cells

    # ---- structures: beta'' first, then the fixtures ----
    structures, summary = [], []
    beta_path = os.path.join(CRYSTALS, "beta_double_prime_Mg5Si6.cif")
    assert os.path.isfile(beta_path), f"run tools/crystal-structures/make_beta_double_prime.py first: {beta_path}"
    e, s = fixture_entry("beta_double_prime", beta_path, "beta_double_prime_Mg5Si6.cif", True, EXTRA_K_MAX)
    structures.append(e); summary.append(s)
    e, s = beta_p1_entry(beta_path)
    structures.append(e); summary.append(s)
    for name in FIXTURE_ORDER:
        path = os.path.join(FIXTURES, name + ".cif")
        assert os.path.isfile(path), path
        e, s = fixture_entry(name, path, name + ".cif", False)
        structures.append(e); summary.append(s)
    out["structures"] = structures

    # ---- the two beta'' entries must agree (detected symmetry equals the file's) ----
    b0, b1 = structures[0], structures[1]
    check("beta'': the P1 list detects the same space group, Laue class, operation count", (
        b0["space_group_number"] == b1["space_group_number"] and b0["laue_class"] == b1["laue_class"]
        and len(b0["operations"]) == len(b1["operations"])))
    check("beta'': the P1 list gives the same allowed count and the same families", (
        b0["reflections"]["allowed_count"] == b1["reflections"]["allowed_count"]
        and [f["name"] for f in b0["families"]] == [f["name"] for f in b1["families"]]))

    # ---- write ----
    os.makedirs(FIXTURES, exist_ok=True)
    dump = lambda o: json.dumps(o, separators=(",", ":"), allow_nan=False)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write('{"versions":' + dump(out["versions"]) + ',\n')
        f.write('"generated_at":' + dump(out["generated_at"]) + ',\n')
        f.write('"definitions":' + dump(out["definitions"]) + ',\n')
        f.write('"cells":[\n' + ",\n".join(dump(c) for c in out["cells"]) + '],\n')
        f.write('"structures":[\n' + ",\n".join(dump(s) for s in out["structures"]) + ']}\n')
    size = os.path.getsize(OUT)
    print(f"\nwrote {OUT}  ({size / 1e6:.2f} MB, {len(cells)} cells, {len(structures)} structures)")

    # ---- registered predictions against what the oracle gave ----
    b = structures[0]
    print(f"\nregistration 5.3 predictions for beta'': allowed 1448 -> {b['reflections']['allowed_count']}, "
          f"families 419 -> {len(b['families'])}; total reflections {b['reflections']['count']}")
    print(f"{len(CHECKS)} self-checks ran; failed: {[n for n, ok in CHECKS if not ok]}")


if __name__ == "__main__":
    main()
