"""DiffractoKit slice 2: the pin generator for scattering factors and structure factors.

Run, from any directory, with the pinned oracle:

    ../References/python/pymatgen/bin/python -I -B tools/pins/structure-factors/pins.py

It writes Tests/DiffractoScatteringTests/Fixtures/numpy-structure-factor-pins.json BEFORE any Swift exists. Deterministic: no clock, no
random numbers, sorted sets; two runs give the same bytes.

Oracle 1 is NumPy, written from the definitions in the comments below (not from any Swift): the table is
np.loadtxt of the fetched py4DSTEM file, F(hkl) = (1/V) sum_sites f_e(Z, |g|^2) occ exp(2 pi i (h x + k y + l z)).
The Laue group (for orbits and names) comes from pymatgen's SpacegroupAnalyzer (symprec 1e-4), the family rule and the
canonical form are slice 1's (copied from tools/pins/cells/pins.py, imported
from nowhere). Oracle 2 is pymatgen's XRDCalculator (an independent implementation of the phase sum with X-ray form
factors). Agreement is agreement, not truth: section 5's closed forms and counts are asserted below too.

The sites the Swift reads
-------------------------
A structure built from a CIF is pinned on the sites AS WRITTEN (pymatgen's CifParser with frac_tolerance = 0, the
expansion slice 1 pinned as `sites_as_written`): the kit's CIFImport keeps the written value (0.333333), pymatgen's
default snaps a coordinate within 1e-4 of 1/3 or 2/3 to the exact fraction. The snapped sites are pinned beside them and
the effect of the snap on F is recorded (`snap_effect`). A presets' sites are exact doubles (1.0/3.0 ...), as the apps'.

Exit status: 0 all self-checks and all predicted counts hold; 1 a self-check failed (assert, loud); 3 the file was
written but a registered PREDICTED COUNT differs from the measured one (printed twice, never "fixed" here).
"""
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
from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS, XRDCalculator
from pymatgen.core import Element, Lattice, Structure
from pymatgen.core.operations import SymmOp
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# ----------------------------------------------------------------------------------------------------------------
# Paths (relative to this script, so it runs from any directory)
# ----------------------------------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
TABLE_PATH = os.path.join(KIT, "References", "py4dstem-scattering", "scattering_factors.txt")
CORE_FIXTURES = os.path.join(KIT, "Tests", "DiffractoCoreTests", "Fixtures")
CRYSTALS = os.path.join(KIT, "References", "crystal-structures")
OUT_DIR = os.path.join(KIT, "Tests", "DiffractoScatteringTests", "Fixtures")
OUT = os.path.join(OUT_DIR, "numpy-structure-factor-pins.json")

# ----------------------------------------------------------------------------------------------------------------
# The registered constants
# ----------------------------------------------------------------------------------------------------------------
TABLE_SHA256 = "be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4"
K_MAX = 1.2                  # Angstrom^-1
EXTRA_K_MAX = (1.6, 2.0)     # beta'' only: kept and ring counts
CUT = 1e-4                   # |F| cut, Angstrom (the apps')
LISTED = 1e-6                # everything above this is listed so the edge band of the cut is visible
LIST_RADIUS = K_MAX * (1 + 1e-6)
FAMILY_TOL = 1e-9            # slice 1's ring rule
TIE_TOL = 1e-12              # slice 1 3.8: equal |g| within 1e-12 relative are ordered by descending (h, k, l)
SYMPREC = 1e-4
ABSENCE_TOL = 1e-9           # slice 1 3.8
WAVELENGTH = 0.5             # Angstrom, the X-ray pin
FE_G2 = (0.0, 0.25, 1.0, 4.0)
HEX_AXES_CLASSES = ("-3", "-3m", "6/m", "6/mmm")

FD_ZS = (1, 6, 12, 13, 14, 29, 74, 79, 103)    # f_e positive and decreasing on [0, 1.2]

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

CHECKS = []          # (name, ok); every check is also an assert
PREDICTIONS = []     # (name, predicted, measured); never asserted, reported at the end


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    assert ok, f"SELF-CHECK FAILED: {name} {detail}"


def predict(name, predicted, measured):
    PREDICTIONS.append((name, predicted, measured))
    if predicted != measured:
        print(f"!!! PREDICTION DIFFERS: {name}: registered {predicted}, measured {measured}", flush=True)


# ----------------------------------------------------------------------------------------------------------------
# The table and f_e
# ----------------------------------------------------------------------------------------------------------------
def load_table():
    data = open(TABLE_PATH, "rb").read()
    got = hashlib.sha256(data).hexdigest()
    check("the table's sha256 is the lock", got == TABLE_SHA256, got)
    t = np.loadtxt(TABLE_PATH, dtype=np.float64)
    check("the table is 103 x 10", t.shape == (103, 10), str(t.shape))
    return t, got


def f_e(table, z, g2):
    """Lobato and Van Dyck, py4DSTEM's units 'A': sum_{i<5} a_i (2 + b_i g2) / (1 + b_i g2)^2, row Z-1, a = columns
    0:10:2, b = columns 1:10:2, g2 = |g|^2 in 1/A^2."""
    row = table[z - 1]
    a, b = row[0:10:2], row[1:10:2]
    g2 = np.asarray(g2, dtype=np.float64)
    return sum(a[i] * (2 + b[i] * g2) / (1 + b[i] * g2) ** 2 for i in range(5))


# ----------------------------------------------------------------------------------------------------------------
# Lattice (slice 1's; the apps' Crystal.latReal), written from the definitions
# ----------------------------------------------------------------------------------------------------------------
def lattice(a, b, c, al, be, ga):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    sg = math.sin(math.radians(ga))
    vol_formula = a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
    A = np.array([[a, 0, 0],
                  [b * cg, b * sg, 0],
                  [c * cb, -c * (cb * cg - ca) / sg, vol_formula / (a * b * sg)]])
    # |det A| by the explicit cofactor expansion (the apps' Crystal.determinant), not LU: exact where the cell is exact
    vol = abs(A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
              - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
              + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
    check(f"cell({a},{b},{c},{al},{be},{ga}): |det A| equals the closed-form volume to 1e-12 relative",
          abs(vol - vol_formula) <= 1e-12 * vol_formula, f"{vol} {vol_formula}")
    L = Lattice.from_parameters(a, b, c, al, be, ga, vesta=True)
    check(f"cell({a},{b},{c},{al},{be},{ga}): A equals pymatgen's vesta matrix to 1e-12",
          np.allclose(A, np.array(L.matrix), rtol=0, atol=1e-12 * max(a, b, c)))
    check(f"cell({a},{b},{c},{al},{be},{ga}): V equals pymatgen's volume to 1e-12 relative",
          abs(vol - L.volume) <= 1e-12 * vol)
    B = np.linalg.inv(A).T          # a_i . b_j = delta_ij; g = h b1 + k b2 + l b3, no 2 pi
    check("a_i . b_j = delta_ij", np.allclose(A @ B.T, np.eye(3), rtol=0, atol=1e-12))
    return A, B, float(vol), L


def metric_is_hexagonal(cell, eps=1e-9):
    a, b, c, al, be, ga = cell
    return abs(a - b) <= eps * a and abs(al - 90) <= eps and abs(be - 90) <= eps and abs(ga - 120) <= eps


# ----------------------------------------------------------------------------------------------------------------
# Slice 1's reflection order, orbits, names, ring rule (copied from tools/pins/cells/pins.py)
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
    out = []
    for x in t:
        k24 = x * 24.0
        k = round(k24)
        check("translation is a multiple of 1/24 within 1e-9", abs(k24 - k) <= 24e-9, f"{x}")
        out.append((k % 24) / 24.0)
    return out


def operations_of(structure):
    """Rotations (int 3x3) and snapped translations as pymatgen's SpacegroupAnalyzer (symprec 1e-4) reports them."""
    sga = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    ops = sga.get_symmetry_operations(cartesian=False)
    rots, trans, xyz = [], [], []
    for op in ops:
        R = np.array(op.rotation_matrix)
        check("rotation part is an integer matrix", np.max(np.abs(R - np.round(R))) < 1e-8)
        R = np.round(R).astype(np.int64)
        t = snap_translation(np.array(op.translation_vector, dtype=float))
        rots.append(R)
        trans.append(t)
        xyz.append(SymmOp.from_rotation_and_translation(R, t).as_xyz_str().replace(" ", ""))
    return sga, np.array(rots), np.array(trans), xyz


def laue_rotations(Rs):
    s = {tuple(int(x) for x in R.ravel()) for R in Rs} | {tuple(int(-x) for x in R.ravel()) for R in Rs}
    return np.array(sorted(s), dtype=np.int64).reshape(-1, 3, 3)


def check_group_closed(name, Rs):
    s = {tuple(int(x) for x in R.ravel()) for R in Rs}
    ok = all(tuple(int(x) for x in (A @ B).ravel()) in s for A in Rs for B in Rs)
    check(f"{name}: Laue rotation set is a group (closed)", ok)


def laue_class_from_group(Rs):
    """Slice 1's Laue class rule, with 'a 4-fold / a 6-fold' read as a proper rotation (det +1) of that order."""
    n = len(Rs)
    dets = [int(round(np.linalg.det(R))) for R in Rs]
    orders = [matrix_order(R) for R in Rs]
    sel = [o for o, d in zip(orders, dets) if d == 1]
    has4, has6 = 4 in sel, 6 in sel
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
    """Slice 1, 3.8: absent iff some operation has (hkl) R == (hkl) and (hkl) . t not an integer within 1e-9."""
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
    """seq: [(hkl tuple, g, |F|^2)] in sweep order (kept reflections). Slice 1's ring rule, plus sumF2."""
    fams = []
    cur = []
    g0 = None

    def close():
        members = sorted((m for m, _, _ in cur), reverse=True)
        orbits = sorted({cache.canonical(m) for m in members}, reverse=True)
        g = min(gg for _, gg, _ in cur)
        fams.append({
            "g": float(g),
            "members": [list(m) for m in members],
            "multiplicity": len(members),
            "orbits": [list(o) for o in orbits],
            "name": " + ".join(orbit_name(o, hex_axes) for o in orbits),
            "sumF2": float(math.fsum(f2 for _, _, f2 in cur)),
        })

    for hkl, g, f2 in seq:
        if cur and g > g0 * (1 + FAMILY_TOL):
            close()
            cur = []
        if not cur:
            g0 = g
        cur.append((hkl, g, f2))
    if cur:
        close()
    return fams


# ----------------------------------------------------------------------------------------------------------------
# Reflections and structure factors, from the definitions
# ----------------------------------------------------------------------------------------------------------------
def enumerate_box(A, B, kmax):
    """Every (h, k, l) != 0 with |g| <= kmax in the per-axis box |h| <= ceil(kmax |a_i|) + 1 (one more than the registered
    bound, so the filter on |g| alone decides membership); returns H (M, 3 int64) and |g| (M,)."""
    n = [int(math.ceil(kmax * np.linalg.norm(A[i]))) + 1 for i in range(3)]
    axes = [np.arange(-m, m + 1) for m in n]
    H = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3).astype(np.int64)
    H = H[np.any(H != 0, axis=1)]
    g = np.linalg.norm(H.astype(float) @ B, axis=1)
    keep = g <= kmax
    return H[keep], g[keep]


def structure_factors(table, vol, H, glen, sites, z_shift=0, power=2, inv_v=True, occ_on=True, form_factor=None):
    """F(hkl) = (1/V) sum_sites f(Z, |g|^2) occ exp(2 pi i (h x + k y + l z)); sites = [(Z, (x, y, z), occ)]."""
    Hf = H.astype(float)
    F = np.zeros(len(H), dtype=np.complex128)
    for (z, pos, occ) in sites:
        fe = form_factor(z, glen ** power) if form_factor else f_e(table, z + z_shift, glen ** power)
        F += fe * (occ if occ_on else 1.0) * np.exp(2j * np.pi * (Hf @ np.array(pos, dtype=float)))
    if inv_v:
        F /= vol
    return F


# ----------------------------------------------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------------------------------------------
FCC = [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)]
_SHIFT = [(x + .25, y + .25, z + .25) for (x, y, z) in FCC]
DIAMOND = FCC + _SHIFT           # the apps' order: the four fcc sites, then the four shifted
HCP = [(0, 0, 0), (2.0 / 3.0, 1.0 / 3.0, 0.5)]
ZS = 0.6225


def sites_of(z, positions, occ=1.0):
    return [(z, tuple(float(x) for x in p), occ) for p in positions]


def preset_structures(table):
    """The apps' presets (Crystal.swift 7d6fdbd1) in the apps' site order, plus macTEM a58c84c's two disordered cases and
    the ideal hcp. Returns [(name, cell, sites, note)]."""
    a_hcp = 3.2094
    q = weak_cut_occupancy(table)
    out = [
        ("Au", (4.0782,) * 3 + (90, 90, 90), sites_of(79, FCC), "apps preset gold: fcc 4.0782"),
        ("Al", (4.0495,) * 3 + (90, 90, 90), sites_of(13, FCC), "apps preset aluminum: fcc 4.0495"),
        ("Si", (5.4309,) * 3 + (90, 90, 90), sites_of(14, DIAMOND), "apps preset silicon: diamond 5.4309, fcc sites then shifted"),
        ("Mg", (3.2094, 3.2094, 5.2108, 90, 90, 120), sites_of(12, HCP), "apps preset magnesium: hcp 3.2094 / 5.2108"),
        ("WS2", (3.1532, 3.1532, 12.323, 90, 90, 120), [
            (74, (1.0 / 3.0, 2.0 / 3.0, 0.25), 1.0),
            (74, (2.0 / 3.0, 1.0 / 3.0, 0.75), 1.0),
            (16, (1.0 / 3.0, 2.0 / 3.0, ZS), 1.0),
            (16, (1.0 / 3.0, 2.0 / 3.0, 1.5 - ZS), 1.0),
            (16, (2.0 / 3.0, 1.0 / 3.0, 1.0 - ZS), 1.0),
            (16, (2.0 / 3.0, 1.0 / 3.0, ZS - 0.5), 1.0),
        ], "apps preset tungstenDisulfide (Schutte 1987), the apps' site order and expressions"),
        ("Cu3Au-occ", (3.75,) * 3 + (90, 90, 90),
         sites_of(79, [(0, 0, 0)], 0.6) + sites_of(29, [(0, 0, 0)], 0.4) + sites_of(29, FCC[1:]),
         "macTEM a58c84c reference.py: Au 0.6 / Cu 0.4 on the origin, Cu on the faces"),
        ("Weak-cut", (3.0,) * 3 + (90, 90, 90), sites_of(29, [(0, 0, 0)]) + sites_of(29, [(.5, .5, .5)], q),
         f"macTEM a58c84c reference.py: CsCl-like Cu, a = 3, second site occupancy q = {q!r} solved so |F(100)| = 2.5e-4"),
        ("GaAs-zb", (5.6533,) * 3 + (90, 90, 90), sites_of(31, FCC) + sites_of(33, _SHIFT),
         "EXTRA (no apps preset): zincblende GaAs, a = 5.6533, Ga on fcc, As on fcc + (1/4, 1/4, 1/4); non-centrosymmetric about the origin: Im F is not zero, its sign is pinned by the closed form"),
        ("ideal_hcp", (a_hcp, a_hcp, a_hcp * math.sqrt(8 / 3), 90, 90, 120), sites_of(12, HCP),
         "ideal hcp, Mg, a = 3.2094, c = a sqrt(8/3)"),
    ]
    # the apps' other presets are fcc/bcc with one element; not pinned here (the apps dump covers them)
    return out, q


def weak_cut_occupancy(table):
    """macTEM a58c84c: q solved so that the lowest odd-sum reflection (1 0 0), |g| = 1/3, has |F| = 2.5e-4:
    F = (1 - q) f_e(29, g^2) / V, V = 27."""
    a = 3.0
    return float(1 - 2.5e-4 * a ** 3 / float(f_e(table, 29, (1 / a) ** 2)))


def raw_cell_from_cif(parser):
    blocks = parser.as_dict()
    assert len(blocks) == 1, "one data block expected"
    data = next(iter(blocks.values()))

    def f(tag):
        return float(re.sub(r"\(.*\)", "", str(data[tag])))

    return tuple(f(t) for t in ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                                "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))


def wrap_unit(x):
    v = x % 1.0
    if v >= 1.0 - 1e-12:
        v = 0.0
    return float(v)


def sites_of_pmg(structure):
    out = []
    for site in structure:
        assert site.is_ordered, "disordered site in a CIF fixture"
        (el, amt), = site.species.items()
        out.append((int(el.Z), tuple(wrap_unit(x) for x in site.frac_coords), float(amt)))
    return out


def parse_cif(path, frac_tolerance=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = CifParser(path) if frac_tolerance is None else CifParser(path, frac_tolerance=frac_tolerance)
        structs = parser.parse_structures(primitive=False)
    assert len(structs) == 1, f"{path}: {len(structs)} structures"
    return parser, structs[0]


CIF_FIXTURES = ["gold_P1", "silicon_Fd-3m", "magnesium_P63mmc", "rocksalt_Fm-3m_192ops", "alpha_quartz_P3221",
                "calcite_R-3c_hex", "fluorapatite_P63m", "pyrite_Pa-3", "ws2_2H", "cdi2_P-3m1"]
EXPECTED_SITES = {"beta_double_prime": 22, "gold_P1": 4, "silicon_Fd-3m": 8, "magnesium_P63mmc": 2,
                  "rocksalt_Fm-3m_192ops": 8, "alpha_quartz_P3221": 9, "calcite_R-3c_hex": 30, "fluorapatite_P63m": 42,
                  "pyrite_Pa-3": 12, "ws2_2H": 6, "cdi2_P-3m1": 3}


def cif_structures():
    """[(name, cell, written sites, snapped sites, symmetry structure (snapped), as-written structure, cif info)]"""
    out = []
    items = [("beta_double_prime", os.path.join(CRYSTALS, "beta_double_prime_Mg5Si6.cif"), "beta_double_prime_Mg5Si6.cif")]
    items += [(n, os.path.join(CORE_FIXTURES, n + ".cif"), n + ".cif") for n in CIF_FIXTURES]
    for name, path, fname in items:
        assert os.path.isfile(path), f"{path} (beta'': run tools/crystal-structures/make_beta_double_prime.py first)"
        text = open(path, encoding="utf-8").read()
        parser, st = parse_cif(path)                       # pymatgen's default: 1/3 and 2/3 snapped
        _, st_w = parse_cif(path, frac_tolerance=0.0)      # as written
        cell = raw_cell_from_cif(parser)
        for x, y in zip(cell, st.lattice.parameters):
            check(f"{name}: cell parameter {x} equals pymatgen's {y} to 1e-12 relative", abs(x - y) <= 1e-12 * max(abs(x), 1.0))
        sn = sites_of_pmg(st)
        wr = sites_of_pmg(st_w)
        check(f"{name}: {EXPECTED_SITES[name]} sites (the cell content the file's formula states), snapped and as written",
              len(sn) == EXPECTED_SITES[name] == len(wr), f"{len(sn)} {len(wr)}")
        out.append((name, cell, wr, sn, st, st_w,
                    {"file": fname, "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(), "text": text}))
    return out


def pmg_structure_from_sites(cell, sites):
    L = Lattice.from_parameters(*cell, vesta=True)
    species = []
    for (z, pos, occ) in sites:
        species.append(Element.from_Z(z).symbol if occ == 1.0 else {Element.from_Z(z).symbol: occ})
    return Structure(L, species, [list(pos) for (z, pos, occ) in sites], coords_are_cartesian=False)


def pmg_structure_cu3au(cell, sites):
    """The two co-sited rows (Au 0.6, Cu 0.4) are ONE disordered site for pymatgen."""
    L = Lattice.from_parameters(*cell, vesta=True)
    species = [{"Au": 0.6, "Cu": 0.4}] + ["Cu"] * 3
    coords = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]
    return Structure(L, species, coords)


# ----------------------------------------------------------------------------------------------------------------
# One structure: reflections, kept, families
# ----------------------------------------------------------------------------------------------------------------
def analyze(S, kmax, detail):
    """S: dict with A, B, vol, sites, table, rot/trans (operations), cache, hex_axes. Returns the listing and families.
    Reflections: every |g| <= kmax (+1e-6 relative when detail) in slice 1's order."""
    radius = LIST_RADIUS if (detail and kmax == K_MAX) else kmax
    H, g = enumerate_box(S["A"], S["B"], radius)
    order = order_reflections(H, g)
    H, g = H[order], g[order]
    F = structure_factors(S["table"], S["vol"], H, g, S["sites"])
    aF = np.abs(F)
    res = {"H": H, "g": g, "F": F, "aF": aF}
    kept = (aF > CUT) & (g <= kmax)
    res["kept"] = kept
    seq = [(tuple(int(x) for x in H[i]), float(g[i]), float(F.real[i] ** 2 + F.imag[i] ** 2)) for i in np.where(kept)[0]]
    res["families"] = build_families(seq, S["cache"], S["hex_axes"])
    allowed = allowed_flags(H, S["rots"], S["trans"])
    res["allowed"] = allowed
    seq_all = [(tuple(int(x) for x in H[i]), float(g[i]), 0.0) for i in range(len(H)) if g[i] <= kmax]
    seq_allowed = [(tuple(int(x) for x in H[i]), float(g[i]), 0.0) for i in range(len(H)) if g[i] <= kmax and allowed[i]]
    res["families_all_n"] = len(build_families(seq_all, S["cache"], S["hex_axes"]))          # every reflection: the lattice rings
    res["families_allowed_n"] = len(build_families(seq_allowed, S["cache"], S["hex_axes"]))  # slice 1: rings of the operation-allowed
    res["allowed_count"] = int((allowed & (g <= kmax)).sum())
    return res


def make_entry(name, how, cell, sites, snapped_sites, sym_structure, build_extra, table, extra_kmax=()):
    A, B, vol, L = lattice(*cell)
    sga, Rs, ts, xyz = operations_of(sym_structure)
    pg = sga.get_point_group_symbol()
    sg_number = int(sga.get_space_group_number())
    check(f"{name}: point group {pg!r} is in the 32 -> 11 table", pg in LAUE_OF_POINT_GROUP)
    LR = laue_rotations(Rs)
    check_group_closed(name, LR)
    laue = laue_class_from_group(LR)
    check(f"{name}: Laue class from the rotation set ({laue}) equals the point group's ({LAUE_OF_POINT_GROUP[pg]})",
          laue == LAUE_OF_POINT_GROUP[pg])
    metric_hex = metric_is_hexagonal(cell)
    hex_axes = laue in HEX_AXES_CLASSES and metric_hex
    S = {"A": A, "B": B, "vol": vol, "sites": sites, "table": table, "rots": Rs, "trans": ts,
         "cache": OrbitCache(LR), "hex_axes": hex_axes, "name": name, "cell": cell, "laue": laue}
    main = analyze(S, K_MAX, True)
    S["main"] = main
    H, g, F, aF = main["H"], main["g"], main["F"], main["aF"]

    # pymatgen's lattice agrees on |g| and the hkl set
    py = L.reciprocal_lattice_crystallographic.get_points_in_sphere([[0.0, 0.0, 0.0]], [0.0, 0.0, 0.0], LIST_RADIUS)
    pys = {tuple(int(round(x)) for x in p[0]) for p in py} - {(0, 0, 0)}
    mine = {tuple(int(x) for x in row) for row in H}
    check(f"{name}: the (h, k, l) set with |g| <= {LIST_RADIUS} equals pymatgen's get_points_in_sphere", pys == mine,
          f"{sorted(pys - mine)[:4]} {sorted(mine - pys)[:4]}")

    # the listing: |F| > LISTED
    listed = aF > LISTED
    items = [[int(H[i][0]), int(H[i][1]), int(H[i][2]), float(g[i]), float(F.real[i]), float(F.imag[i]), float(aF[i])]
             for i in np.where(listed)[0]]
    kept = main["kept"]
    kept_count = int(kept.sum())

    # operation-absent reflections must have |F| <= 1e-12 (on the exact sites: preset, or the snapped CIF sites)
    allowed = main["allowed"]
    absent = ~allowed
    exact_sites = snapped_sites if snapped_sites is not None else sites
    F_exact = structure_factors(table, vol, H, g, exact_sites)
    absent_max_exact = float(np.abs(F_exact[absent]).max()) if absent.any() else 0.0
    absent_max_written = float(aF[absent].max()) if absent.any() else 0.0
    check(f"{name}: every operation-absent reflection ({int(absent.sum())}) has |F| <= 1e-12 (exact sites)",
          absent_max_exact <= 1e-12, f"{absent_max_exact}")
    if snapped_sites is not None:
        snap_effect = float(np.max(np.abs(np.abs(F_exact) - aF) / np.maximum(aF, CUT)))
        kept_snapped = (np.abs(F_exact) > CUT) & (g <= K_MAX)
        snap_kept_differs = int((kept_snapped != ((aF > CUT) & (g <= K_MAX))).sum())
    else:
        snap_effect = None
        snap_kept_differs = None
    # a single orbit is wholly kept or wholly dropped
    cache = S["cache"]
    seen = {}
    for i in range(len(H)):
        key = cache.canonical(tuple(int(x) for x in H[i]))
        seen.setdefault(key, []).append(i)
    split = []
    for key, idx in seen.items():
        within = [i for i in idx if g[i] <= K_MAX]
        if not within:
            continue
        flags = {bool(kept[i]) for i in within}
        if len(flags) > 1:
            near = [i for i in within if abs(aF[i] - CUT) <= 1e-12 * CUT]
            if near:
                continue
            split.append((key, [(tuple(int(x) for x in H[i]), float(aF[i])) for i in within]))
    check(f"{name}: a single orbit is wholly kept or wholly dropped", not split, f"{split[:2]}")

    entry = {
        "name": name,
        "build": dict({"how": how}, **build_extra),
        "cell": {"a": cell[0], "b": cell[1], "c": cell[2], "alpha": cell[3], "beta": cell[4], "gamma": cell[5]},
        "volume": vol,
        "sites": [{"z": int(z), "frac": [float(x) for x in pos], "occupancy": float(occ)} for (z, pos, occ) in sites],
        "space_group_number": sg_number,
        "point_group": pg,
        "laue_class": laue,
        "metric_hexagonal": bool(metric_hex),
        "uses_miller_bravais": bool(hex_axes),
        "laue_rotations": LR.tolist(),
        "operation_count": int(len(Rs)),
        "reflections": {"k_max": K_MAX, "list_radius": LIST_RADIUS, "listed_floor": LISTED,
                        "total_in_radius": int(len(H)), "listed_count": len(items), "items": items},
        "kept_count": kept_count,
        "families": main["families"],
        "ring_counts": {"after_cut": len(main["families"]), "operation_allowed": main["families_allowed_n"],
                        "all_reflections": main["families_all_n"]},
        "operation_allowed_count": main["allowed_count"],
        "max_ImF_over_absF_kept": float(np.max(np.abs(F.imag[kept]) / aF[kept])) if kept.any() else 0.0,
        "operation_absent": {"count": int(absent.sum()), "max_absF_exact_sites": absent_max_exact,
                             "max_absF_pinned_sites": absent_max_written},
        "cut_edge_band": {"kept_in_1e-4_1e-3": int(((aF > CUT) & (aF < 1e-3) & (g <= K_MAX)).sum()),
                          "dropped_in_1e-5_1e-4": int(((aF <= CUT) & (aF > 1e-5) & (g <= K_MAX)).sum()),
                          "closest_to_cut_relative": float(np.min(np.abs(aF[g <= K_MAX] - CUT) / CUT))
                          if (g <= K_MAX).any() else None},
    }
    if snapped_sites is not None:
        entry["sites_snapped"] = [{"z": int(z), "frac": [float(x) for x in pos], "occupancy": float(occ)}
                                  for (z, pos, occ) in snapped_sites]
        entry["snap_effect"] = {"max_relative_absF_difference": snap_effect, "kept_set_differences": snap_kept_differs,
                                "note": "sites as written (pinned) vs pymatgen's 1/3, 2/3 snapped sites"}
    extras = {}
    for km in extra_kmax:
        r = analyze(S, km, False)
        extras[repr(km)] = {"total": int(len(r["H"])), "kept_count": int(r["kept"].sum()),
                            "rings": len(r["families"]), "rings_operation_allowed": r["families_allowed_n"],
                            "rings_all_reflections": r["families_all_n"], "operation_allowed_count": r["allowed_count"]}
    if extras:
        entry["extra_k_max"] = extras
    S["entry"] = entry
    print(f"\n=== {name}  [{how}]  sg {sg_number} {pg}  Laue {laue}{' (Miller-Bravais names)' if hex_axes else ''}  "
          f"sites {len(sites)}  V = {vol:.6f}")
    print(f"  reflections in radius {len(H)}, listed (|F| > 1e-6) {len(items)}, kept {kept_count}, rings after the cut "
          f"{len(main['families'])} (operation-allowed {main['families_allowed_n']}, all {main['families_all_n']}); operation-absent {int(absent.sum())} "
          f"(max |F| exact {absent_max_exact:.2e}, pinned {absent_max_written:.2e})"
          + (f"; snap effect {snap_effect:.2e} (kept sets differ in {snap_kept_differs})" if snap_effect is not None else ""))
    for km, r in extras.items():
        print(f"  at {km}: kept {r['kept_count']}, rings {r['rings']} (operation-allowed {r['rings_operation_allowed']}, all {r['rings_all_reflections']})")
    return S


# ----------------------------------------------------------------------------------------------------------------
# Truth: extinctions, closed forms, weak cut
# ----------------------------------------------------------------------------------------------------------------
def truth_checks(table, SS):
    def S_of(n):
        return SS[n]

    def parts(S):
        m = S["main"]
        return m["H"], m["g"], m["F"], m["aF"]

    # f_e itself
    gs = np.linspace(0, K_MAX, 241) ** 2
    for z in FD_ZS:
        fe = f_e(table, z, gs)
        check(f"f_e(Z={z}) is positive and strictly decreasing on [0, 1.2]", bool((fe > 0).all() and (np.diff(fe) < 0).all()))
        row = table[z - 1]
        check(f"f_e(Z={z}, 0) = 2 sum a_i", abs(float(f_e(table, z, 0.0)) - 2 * float(np.sum(row[0:10:2]))) <= 1e-13 * abs(2 * float(np.sum(row[0:10:2]))) + 1e-15)

    def rel_close(computed, closed, scale, label):
        """1e-13 relative to |closed|; where |closed| < 1e-3 of its scale, 1e-13 of the scale (absolute)."""
        err = np.abs(computed - closed)
        tol = 1e-13 * np.where(np.abs(closed) > 1e-3 * scale, np.abs(closed), scale)
        check(label, bool((err <= tol).all()), f"worst {float((err / tol).max())} of the tolerance")
        return float(np.max(err / np.maximum(np.abs(closed), 1e-3 * scale)))

    out = {}
    # fcc: Au, Al
    for n in ("Au", "Al"):
        S = S_of(n)
        H, g, F, aF = parts(S)
        mixed = np.ptp(H % 2, axis=1) != 0
        zz = S["sites"][0][0]
        fe = f_e(table, zz, g ** 2)
        check(f"{n}: |F| <= 1e-12 for h, k, l of mixed parity", bool((aF[mixed] <= 1e-12).all()), f"{aF[mixed].max()}")
        check(f"{n}: |F| > 1e-3 for unmixed within 1.2", bool((aF[~mixed & (g <= K_MAX)] > 1e-3).all()), f"{aF[~mixed].min()}")
        closed = 4 * fe[~mixed] / S["vol"]
        out[f"{n} fcc 4 f/V"] = rel_close(F[~mixed], closed.astype(complex), 4 * fe[~mixed].max() / S["vol"],
                                           f"{n}: F = 4 f_e / V for unmixed hkl, 1e-13 relative")
    # diamond: Si
    S = S_of("Si")
    H, g, F, aF = parts(S)
    fe = f_e(table, 14, g ** 2)
    par = np.ptp(H % 2, axis=1)
    mixed = par != 0
    fourn2 = (~mixed) & (H.sum(axis=1) % 2 == 0) & (H.sum(axis=1) % 4 == 2)
    check("Si: |F| <= 1e-12 for mixed parity", bool((aF[mixed] <= 1e-12).all()), f"{aF[mixed].max()}")
    check("Si: |F| <= 1e-12 for all-even h + k + l = 4n + 2 ({2 0 0}, {2 2 2}, {4 2 0}, {4 4 2}, ...)", bool((aF[fourn2] <= 1e-12).all()),
          f"{aF[fourn2].max()}")
    sample = {tuple(sorted((abs(int(x)) for x in H[i]), reverse=True)) for i in np.where(fourn2 & (g <= K_MAX))[0]}
    check("Si: the extinct {2 0 0} {2 2 2} {4 2 0} {4 4 2} are among the fourn+2 reflections within 1.2",
          {(2, 0, 0), (2, 2, 2), (4, 2, 0), (4, 4, 2)} <= sample, f"{sorted(sample)}")
    allowed_si = (~mixed) & (~fourn2)
    check("Si: |F| > 1e-3 for the others within 1.2", bool((aF[allowed_si & (g <= K_MAX)] > 1e-3).all()),
          f"{aF[allowed_si & (g <= K_MAX)].min()}")
    ssum = H.sum(axis=1)
    i_pow = np.array([1, 1j, -1, -1j])[ssum % 4]
    closed = 4 * fe[allowed_si] * (1 + i_pow[allowed_si]) / S["vol"]
    out["Si diamond 4 f (1 + i^(h+k+l)) / V"] = rel_close(F[allowed_si], closed, 8 * fe[allowed_si].max() / S["vol"],
                                                           "Si: F = 4 f_e (1 + i^(h+k+l)) / V for unmixed hkl, 1e-13 relative")
    odd = allowed_si & (ssum % 2 != 0)
    even = allowed_si & (ssum % 4 == 0)
    check("Si: |F| = 4 sqrt(2) f / V for h + k + l odd and 8 f / V for h + k + l = 4n",
          bool(np.allclose(aF[odd], 4 * math.sqrt(2) * fe[odd] / S["vol"], rtol=1e-13, atol=0)
               and np.allclose(aF[even], 8 * fe[even] / S["vol"], rtol=1e-13, atol=0)))
    # zincblende: GaAs (the Im F sign is pinned by the closed form)
    S = S_of("GaAs-zb")
    H, g, F, aF = parts(S)
    fGa, fAs = f_e(table, 31, g ** 2), f_e(table, 33, g ** 2)
    mixed = np.ptp(H % 2, axis=1) != 0
    check("GaAs-zb: |F| <= 1e-12 for mixed parity", bool((aF[mixed] <= 1e-12).all()), f"{aF[mixed].max()}")
    ssum = H.sum(axis=1)
    i_pow = np.array([1, 1j, -1, -1j])[ssum % 4]
    um = ~mixed
    closed = 4 * (fGa[um] + i_pow[um] * fAs[um]) / S["vol"]
    out["GaAs-zb 4 (fGa + i^(h+k+l) fAs) / V"] = rel_close(F[um], closed, 4 * (fGa[um] + fAs[um]).max() / S["vol"],
                                                         "GaAs-zb: F = 4 (f_Ga + i^(h+k+l) f_As) / V for unmixed hkl, 1e-13 relative (pins the sign of Im F)")
    h111 = (H == [1, 1, 1]).all(axis=1)
    check("GaAs-zb: Im F(1 1 1) < 0 and Im F(-1 -1 -1) > 0 (i^3 = -i), the sign a phase flip would change",
          bool(F[h111][0].imag < 0 and F[(H == [-1, -1, -1]).all(axis=1)][0].imag > 0))
    # hcp: Mg
    S = S_of("Mg")
    H, g, F, aF = parts(S)
    h, k, l = H[:, 0], H[:, 1], H[:, 2]
    absent = (l % 2 == 1) & (((h - k) % 3) == 0)
    present = ~absent
    check("Mg: |F| <= 1e-12 where l is odd and h - k = 3n (includes 000l, l odd)", bool((aF[absent] <= 1e-12).all()), f"{aF[absent].max()}")
    check("Mg: 000l with l odd is among them", bool(((h == 0) & (k == 0) & (l % 2 == 1) & absent).any()))
    check("Mg: |F| > 1e-3 elsewhere within 1.2", bool((aF[present & (g <= K_MAX)] > 1e-3).all()), f"{aF[present & (g <= K_MAX)].min()}")
    fe = f_e(table, 12, g ** 2)
    closed = 2 * fe * np.abs(np.cos(np.pi * ((h + 2 * k) / 3 + l / 2))) / S["vol"]
    err = np.abs(aF - closed)
    scale = 2 * fe / S["vol"]
    check("Mg: |F| = 2 f |cos(pi((h + 2k)/3 + l/2))| / V, 1e-13 relative (of the scale where the closed form is below 1e-3 of it)",
          bool((err <= 1e-13 * np.where(closed > 1e-3 * scale, closed, scale)).all()), f"{float(np.max(err / scale))}")
    out["Mg hcp"] = float(np.max(err / np.maximum(closed, 1e-3 * scale)))
    # Cu3Au-occ
    S = S_of("Cu3Au-occ")
    H, g, F, aF = parts(S)
    fAu, fCu = f_e(table, 79, g ** 2), f_e(table, 29, g ** 2)
    mixed = np.ptp(H % 2, axis=1) != 0
    closed = np.where(mixed, (0.6 * fAu + 0.4 * fCu - fCu), (0.6 * fAu + 3.4 * fCu)) / S["vol"]
    out["Cu3Au-occ"] = rel_close(F, closed.astype(complex), float(np.max(np.abs(closed))),
                                 "Cu3Au-occ: unmixed (0.6 fAu + 3.4 fCu)/V, mixed (0.6 fAu + 0.4 fCu - fCu)/V, 1e-13 relative")
    check("Cu3Au-occ: the superlattice reflections (mixed parity) exist only through the occupancy: |F| > 1e-3 and < the fundamentals",
          bool((aF[mixed & (g <= K_MAX)] > 1e-3).all()))
    # Weak-cut straddle
    S = S_of("Weak-cut")
    H, g, F, aF = parts(S)
    fe = f_e(table, 29, g ** 2)
    q = S["sites"][1][2]
    oddsum = H.sum(axis=1) % 2 != 0
    closed = np.where(oddsum, (1 - q), (1 + q)) * fe / S["vol"]
    out["Weak-cut"] = rel_close(F, closed.astype(complex), float(np.max(np.abs(closed))),
                                "Weak-cut: F = (1 -/+ q) f_e / V for h + k + l odd / even, 1e-13 relative")
    inb = g <= K_MAX
    nk = int(((aF > CUT) & (aF < 1e-3) & inb & oddsum).sum())
    nd = int(((aF <= CUT) & (aF > 1e-5) & inb & oddsum).sum())
    check("Weak-cut: at least 3 reflections kept in (1e-4, 1e-3) and 3 dropped in (1e-5, 1e-4)", nk >= 3 and nd >= 3, f"{nk} {nd}")
    kept_side = bool(((closed > CUT) == (aF > CUT))[inb].all())
    check("Weak-cut: every reflection is on the side of the cut that the closed form puts it", kept_side)
    out["Weak-cut straddle"] = {"kept_in_1e-4_1e-3": nk, "dropped_in_1e-5_1e-4": nd, "q": q,
                                "|F(100)|": float(aF[(H == [1, 0, 0]).all(axis=1)][0])}
    check("Weak-cut: |F(100)| = 2.5e-4 by the solved occupancy to 1e-12 relative", abs(out["Weak-cut straddle"]["|F(100)|"] - 2.5e-4) <= 2.5e-16)
    return out


# ----------------------------------------------------------------------------------------------------------------
# Mutation-moves: each planned mutation must move its target in NumPy by more than 1e-6
# ----------------------------------------------------------------------------------------------------------------
def mutation_moves(table, SS):
    moves = {}      # mutation -> {structure: value}
    designated = {
        "table_row_z_plus1": ["Au"],
        "fe_of_g_not_g2": ["Al"],
        "no_inverse_volume": ["Mg"],
        "occupancy_ignored": ["Cu3Au-occ", "Weak-cut"],
        "phase_sign_flipped_ImF": ["alpha_quartz_P3221", "GaAs-zb"],
        "diamond_basis_dropped": ["Si"],
        "last_site_skipped": ["Au", "Si", "WS2"],
        "cut_1e-3_kept_count": ["Weak-cut"],
        "families_on_operation_allowed_list": ["Si", "Mg"],
    }

    def relabs(F2, S):
        m = S["main"]
        k = m["kept"]
        return float(np.max(np.abs(np.abs(F2[k]) - m["aF"][k]) / m["aF"][k]))

    for name, S in SS.items():
        m = S["main"]
        H, g, aF, F = m["H"], m["g"], m["aF"], m["F"]
        sites, vol = S["sites"], S["vol"]
        kept = m["kept"]

        def rec(mut, val):
            moves.setdefault(mut, {})[name] = float(val)

        zs = [z for (z, _, _) in sites]
        if max(zs) + 1 <= 103:
            rec("table_row_z_plus1", relabs(structure_factors(table, vol, H, g, sites, z_shift=+1), S))
        if min(zs) - 1 >= 1:
            rec("table_row_z_minus1", relabs(structure_factors(table, vol, H, g, sites, z_shift=-1), S))
        rec("fe_of_g_not_g2", relabs(structure_factors(table, vol, H, g, sites, power=1), S))
        rec("no_inverse_volume", relabs(structure_factors(table, vol, H, g, sites, inv_v=False), S))
        rec("occupancy_ignored", relabs(structure_factors(table, vol, H, g, sites, occ_on=False), S))
        # phase sign flipped: the conjugate; |F| does not move, Im F does (relative to |F|)
        Fc = np.conj(F)
        rec("phase_sign_flipped_absF", relabs(Fc, S))
        rec("phase_sign_flipped_ImF", np.max(np.abs(Fc.imag[kept] - F.imag[kept]) / aF[kept]))
        if len(sites) > 1:
            rec("last_site_skipped", relabs(structure_factors(table, vol, H, g, sites[:-1]), S))
        # cut at 1e-3: how many kept reflections change
        rec("cut_1e-3_kept_count", int(kept.sum()) - int(((aF > 1e-3) & (g <= K_MAX)).sum()))
        rec("cut_ge_reflections_exactly_at_cut", int((aF[g <= K_MAX] == CUT).sum()))
        rec("families_on_operation_allowed_list", m["families_allowed_n"] - len(m["families"]))
        rec("families_on_all_reflections", m["families_all_n"] - len(m["families"]))
        if name == "Si":
            sub = sites_of(14, FCC)
            F2 = structure_factors(table, vol, H, g, sub)
            allowed = aF > 1e-3
            rec("diamond_basis_dropped", float(np.max(np.abs(np.abs(F2[allowed]) - aF[allowed]) / aF[allowed])))
            fourn2 = (np.ptp(H % 2, axis=1) == 0) & (H.sum(axis=1) % 2 == 0) & (H.sum(axis=1) % 4 == 2)
            rec("diamond_basis_dropped_4n+2_leak_absF", float(np.abs(F2[fourn2]).max()))

    for mut, names in designated.items():
        for n in names:
            v = moves[mut][n]
            check(f"mutation-move: '{mut}' moves its target on {n} by {v:.3e} (> 1e-6)", v > 1e-6)
    check("mutation-move: the phase-conjugate leaves |F| unmoved on WS2 (< 1e-12): only Im F and the X-ray pin see it",
          moves["phase_sign_flipped_absF"]["WS2"] < 1e-12)
    check("mutation-move: Si's basis-dropped leak into the 4n+2 reflections exceeds 1e-3",
          moves["diamond_basis_dropped_4n+2_leak_absF"]["Si"] > 1e-3)
    return moves, designated


# ----------------------------------------------------------------------------------------------------------------
# pymatgen's X-ray pattern
# ----------------------------------------------------------------------------------------------------------------
def xray_f(z, s2):
    sym = Element.from_Z(z).symbol
    coeff = np.asarray(ATOMIC_SCATTERING_PARAMS[sym])
    return z - 41.78214 * s2 * np.sum(coeff[:, 0] * np.exp(-coeff[:, 1] * s2[:, None]), axis=1)


def xray_block(name, S, pmg_struct):
    calc = XRDCalculator(wavelength=WAVELENGTH)
    # two_theta_range=None (every beam in the limiting sphere, |g| <= 2/lambda = 4) lets the Lorentz-polarisation factor,
    # which diverges as 2theta -> 180, set the maximum, and pymatgen then drops every low-angle peak below its 1e-3 percent
    # floor (Cu3Au-occ: all of them). The range stops 0.1 percent beyond the pinned radius instead.
    tt_max = 2 * math.degrees(math.asin(WAVELENGTH * LIST_RADIUS * (1 + 1e-3) / 2))
    pat = calc.get_pattern(pmg_struct, scaled=False, two_theta_range=(0, tt_max))
    is_hex = pmg_struct.lattice.is_hexagonal()
    A, B = S["A"], S["B"]
    H, g = enumerate_box(A, B, LIST_RADIUS)
    half = (H[:, 0] > 0) | ((H[:, 0] == 0) & (H[:, 1] > 0)) | ((H[:, 0] == 0) & (H[:, 1] == 0) & (H[:, 2] > 0))
    H, g = H[half], g[half]
    order = order_reflections(H, g)
    H, g = H[order], g[order]
    # my X-ray structure factors for every half-space reflection
    s2 = (g / 2) ** 2
    Hf = H.astype(float)
    F = np.zeros(len(H), dtype=np.complex128)
    for (z, pos, occ) in S["sites"]:
        F += xray_f(z, s2) * occ * np.exp(2j * np.pi * (Hf @ np.array(pos)))
    absF2 = F.real ** 2 + F.imag ** 2
    theta_all = np.arcsin(np.clip(WAVELENGTH * g / 2, -1, 1))
    lp_all = (1 + np.cos(2 * theta_all) ** 2) / (np.sin(theta_all) ** 2 * np.cos(theta_all))
    my_int = 2 * absF2 * lp_all
    max_int = max(pat.y)

    peaks = []
    matched = np.zeros(len(H), dtype=bool)
    worst = 0.0
    for x, y, hk, d in zip(pat.x, pat.y, pat.hkls, pat.d_hkls):
        gp = 1.0 / d
        if gp > LIST_RADIUS:
            continue
        theta = math.radians(x) / 2
        lp = (1 + math.cos(2 * theta) ** 2) / (math.sin(theta) ** 2 * math.cos(theta))
        sum_f2_pmg = y / (2 * lp)
        mem = np.where(np.abs(g - gp) <= 1e-9 * gp)[0]
        check(f"{name}: X-ray peak at 2theta {x:.6f} has members", len(mem) > 0)
        matched[mem] = True
        mine = float(math.fsum(absF2[mem]))
        rel = abs(mine - sum_f2_pmg) / sum_f2_pmg
        worst = max(worst, rel)
        peaks.append({"two_theta": float(x), "intensity": float(y), "d": float(d), "g": float(gp),
                      "pymatgen_hkls": [{"hkl": [int(v) for v in e["hkl"]], "multiplicity": int(e["multiplicity"])} for e in hk],
                      "lp": lp, "sum_F2_pymatgen": float(sum_f2_pmg),
                      "members": [[int(v) for v in H[i]] for i in mem],
                      "sum_F2_numpy": mine, "relative_difference": rel})
    check(f"{name}: X-ray sum of |F_X|^2 over the half-space members equals pymatgen's within 1e-6 relative, every peak "
          f"(worst {worst:.2e}, {len(peaks)} peaks)", worst <= 1e-6)
    # every half-space reflection pymatgen did not report has a negligible intensity (below its scaled tolerance)
    unmatched_int = my_int[~matched & (g <= LIST_RADIUS)]
    worst_unmatched = float(unmatched_int.max() / max_int * 100) if len(unmatched_int) else 0.0
    check(f"{name}: reflections absent from pymatgen's peak list are below its 1e-3 percent intensity floor (largest {worst_unmatched:.2e} percent)",
          worst_unmatched <= 1e-3 * (1 + 1e-6))
    elems = sorted({z for (z, _, _) in S["sites"]})
    coeffs = {Element.from_Z(z).symbol: {"z": z, "coefficients": [list(map(float, c)) for c in ATOMIC_SCATTERING_PARAMS[Element.from_Z(z).symbol]]}
              for z in elems}
    print(f"  X-ray {name}: {len(peaks)} peaks within 1.2, worst relative difference {worst:.2e}; "
          f"largest unreported peak {worst_unmatched:.2e} percent of max")
    return {"name": name, "wavelength": WAVELENGTH, "two_theta_range": [0.0, tt_max], "is_hexagonal_lattice": bool(is_hex),
            "hkl_note": "pymatgen_hkls are pymatgen's own labels (four indices h k i l in a hexagonal lattice); members are the (h, k, l) half-space members used for the sums",
            "element_coefficients": coeffs, "peak_count": len(peaks), "worst_relative_difference": worst,
            "largest_unreported_percent_of_max": worst_unmatched, "peaks": peaks}


# ----------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------
def main():
    table, table_sha = load_table()
    versions = {"python": sys.version, "pymatgen": dist_version("pymatgen"), "pymatgen_core": dist_version("pymatgen-core"),
                "spglib": spglib.__version__, "numpy": np.__version__}
    check("pymatgen.core.__version__ equals the pymatgen-core distribution", pymatgen.core.__version__ == versions["pymatgen_core"])
    print("oracle:", json.dumps(versions))

    out = {
        "versions": versions,
        "table_sha256": table_sha,
        "definitions": {
            "cut": CUT, "listed_floor": LISTED, "k_max": K_MAX, "list_radius": LIST_RADIUS, "family_tolerance": FAMILY_TOL,
            "tie_tolerance": TIE_TOL, "symprec": SYMPREC, "xray_wavelength": WAVELENGTH, "extra_k_max": list(EXTRA_K_MAX),
            "fe": "sum_{i<5} a_i (2 + b_i g2) / (1 + b_i g2)^2, a = row[0:10:2], b = row[1:10:2], row Z-1, g2 = |g|^2 in 1/A^2",
            "F": "(1/V) sum_sites f_e(Z, |g|^2) occ exp(2 pi i (h x + k y + l z)), V = |det A|, A rows (a,0,0), (b cos g, b sin g, 0), (c cos b, -c (cos b cos g - cos a)/sin g, V/(a b sin g)), B = inv(A)^T",
            "kept": "|F| > cut (strict) and |g| <= k_max",
            "listed": "|F| > listed_floor and |g| <= list_radius = k_max (1 + 1e-6)",
            "reflection_order": "increasing |g|; ties (within 1e-12 relative of the cluster's first) by descending (h, k, l)",
            "families": "slice 1 rule on the kept reflections: join while |g| <= g0 (1 + 1e-9); g = smallest member; orbits by the Laue group (rotation set +-R from SpacegroupAnalyzer symprec 1e-4), canonical = max by (count of non-negative indices, lexicographic); name = orbits descending joined by ' + ', Miller-Bravais only for a hexagonal-axes Laue class in a hexagonal metric; sumF2 = sum of Re^2 + Im^2 over the members (math.fsum)",
            "cif_sites": "sites AS WRITTEN (pymatgen CifParser frac_tolerance = 0, expanded, wrapped to [0, 1)); snapped sites (default CifParser) pinned beside them as sites_snapped",
        },
    }

    # ---- the table and f_e samples ----
    out["table"] = table.tolist()
    fe_rows = []
    for z in range(1, 104):
        fe_rows.append([float(f_e(table, z, g2)) for g2 in FE_G2])
    out["fe_samples"] = {"g2": list(FE_G2), "rows": fe_rows, "note": "rows[Z-1][i] = f_e(Z, g2[i]) in Angstrom, NumPy"}
    check("f_e samples are finite", bool(np.isfinite(np.array(fe_rows)).all()))

    # ---- structures ----
    presets, q = preset_structures(table)
    SS = {}
    order = []
    for name, cell, sites, note in presets:
        if name == "Cu3Au-occ":
            sym = pmg_structure_cu3au(cell, sites)
        else:
            sym = pmg_structure_from_sites(cell, sites)
        S = make_entry(name, "preset", cell, sites, None, sym, {"note": note}, table)
        S["pmg"] = sym
        SS[name] = S
        order.append(name)
    for name, cell, wr, sn, st_snapped, st_written, info in cif_structures():
        sga_w = SpacegroupAnalyzer(st_written, symprec=SYMPREC)
        sga_s = SpacegroupAnalyzer(st_snapped, symprec=SYMPREC)
        check(f"{name}: the sites as written give the same space group and operation count as the snapped ones",
              sga_w.get_space_group_number() == sga_s.get_space_group_number()
              and len(sga_w.get_symmetry_operations()) == len(sga_s.get_symmetry_operations()))
        extra = EXTRA_K_MAX if name == "beta_double_prime" else ()
        build = {"cif_file": info["file"], "cif_sha256": info["sha256"],
                 "swift": "built through CIFImport from this file's text (cif_text for beta'', the fixture file otherwise)"}
        if name == "beta_double_prime":
            build["cif_text"] = info["text"]
        S = make_entry(name, "cif", cell, wr, sn, st_snapped, build, table, extra)
        S["pmg"] = st_written
        SS[name] = S
        order.append(name)

    # ---- truth ----
    print("\n--- truth: extinctions, closed forms, the weak cut ---")
    truth = truth_checks(table, SS)
    print(json.dumps(truth, default=float, indent=1))

    # ---- registered predicted counts ----
    print("\n--- predicted counts (registration 5, 'rings after the cut') ---")
    def rings(n): return len(SS[n]["entry"]["families"])
    predict("Au rings after the cut", 8, rings("Au"))
    predict("Al rings after the cut", 8, rings("Al"))
    predict("Si rings after the cut", 10, rings("Si"))
    def rc(n, k): return SS[n]["entry"]["ring_counts"][k]
    predict("Si rings by symmetry (operation-allowed)", 12, rc("Si", "operation_allowed"))
    predict("Mg rings after the cut", 24, rings("Mg"))
    predict("Mg rings by symmetry (operation-allowed)", 25, rc("Mg", "operation_allowed"))
    predict("ideal hcp rings after the cut", 23, rings("ideal_hcp"))
    predict("ideal hcp rings by symmetry (operation-allowed)", 24, rc("ideal_hcp", "operation_allowed"))
    b = SS["beta_double_prime"]["entry"]
    predict("beta'' rings after the cut at 1.2", 417, len(b["families"]))
    predict("beta'' kept reflections at 1.2", 1442, b["kept_count"])
    predict("beta'' rings after the cut at 1.6", 953, b["extra_k_max"]["1.6"]["rings"])
    predict("beta'' rings after the cut at 2.0", 1805, b["extra_k_max"]["2.0"]["rings"])
    predict("beta'' operation-allowed reflections at 1.2 (slice 1)", 1448, b["operation_allowed_count"])
    predict("beta'' operation-allowed rings at 1.2 (slice 1)", 419, rc("beta_double_prime", "operation_allowed"))
    members_b = sum(f["multiplicity"] for f in b["families"])
    predict("beta'' reflections in the rings at 1.2 (sum of multiplicities)", 1442, members_b)
    # operation-allowed reflections the cut dropped in beta'' (the prediction names {3 3 3} and {14 0 1})
    m = SS["beta_double_prime"]["main"]
    drop_names = {}
    for i in range(len(m["H"])):
        if m["allowed"][i] and m["g"][i] <= K_MAX and not m["kept"][i]:
            c = SS["beta_double_prime"]["cache"].canonical(tuple(int(x) for x in m["H"][i]))
            drop_names.setdefault(c, []).append(float(m["aF"][i]))
    print("beta'' operation-allowed reflections dropped by the cut, by orbit:",
          {orbit_name(c, False): (len(v), f"{max(v):.2e}") for c, v in sorted(drop_names.items())})
    predict("beta'' operation-allowed orbits dropped by the cut: {3 3 3} (4, 9.3e-5) and {14 0 1} (2, 7.9e-5)",
            {"{3 3 3}": (4, 9.3), "{14 0 1}": (2, 7.9)},
            {orbit_name(c, False): (len(v), round(max(v) * 1e5, 1)) for c, v in sorted(drop_names.items())})
    out["beta_double_prime_dropped_orbits"] = {orbit_name(c, False): {"members": len(v), "max_absF": max(v)} for c, v in sorted(drop_names.items())}

    # ---- mutation-moves ----
    print("\n--- mutation-moves: each planned mutant must move its target in NumPy by more than 1e-6 ---")
    moves, designated = mutation_moves(table, SS)
    for mut in sorted(moves):
        tag = ", ".join(f"{n} {v:.3e}" for n, v in moves[mut].items() if n in designated.get(mut, []))
        if tag:
            print(f"  {mut}: {tag}")
    out["mutation_moves"] = {"threshold": 1e-6, "designated": designated, "moves": moves,
                             "notes": {
                                 "table_row_z_plus1": "relative change of |F| over the kept reflections, reading row Z instead of row Z-1",
                                 "fe_of_g_not_g2": "relative change of |F|, f_e(|g|) instead of f_e(|g|^2)",
                                 "no_inverse_volume": "relative change of |F| without the 1/V",
                                 "occupancy_ignored": "relative change of |F| with occ = 1 for every site",
                                 "phase_sign_flipped_ImF": "max |ImF' - ImF| / |F| over the kept reflections; |F| itself never moves (phase_sign_flipped_absF), and neither does sum |F|^2 over a ring, so the X-ray pin cannot see it. WS2 (the apps' preset and the CIF) is centrosymmetric about the cell origin (W at +-(1/3, 2/3, 1/4)): Im F = 0 there, so the flip is invisible on it; quartz and GaAs-zb see it",
                                 "diamond_basis_dropped": "Si with the fcc sites only: relative change of |F| on the allowed reflections; the 4n+2 leak is the largest |F| it puts on reflections that must be absent",
                                 "last_site_skipped": "relative change of |F|",
                                 "cut_1e-3_kept_count": "kept at 1e-4 minus kept at 1e-3 (a count)",
                                 "cut_ge_reflections_exactly_at_cut": "number of reflections with |F| == 1e-4 exactly: 0 for every pinned structure, so '>=' for '>' is invisible on them; see cut_edge",
                                 "families_on_operation_allowed_list": "rings of the operation-allowed reflections (slice 1) minus rings after the cut (a count)",
                                 "families_on_all_reflections": "rings of every reflection (the lattice rings) minus rings after the cut (a count)",
                             }}

    # ---- the cut edge: a cell whose reflections all have |F| == 1e-4 exactly (a constant unit form factor) ----
    ce_cell = (100.0, 100.0, 1.0, 90, 90, 90)
    A, B, vol, _ = lattice(*ce_cell)
    Hc, gc = enumerate_box(A, B, K_MAX)
    Fc = structure_factors(table, vol, Hc, gc, [(1, (0.0, 0.0, 0.0), 1.0)], form_factor=lambda z, g2: np.ones_like(g2))
    re_, im_ = Fc.real, Fc.imag
    absF = np.sqrt(re_ * re_ + im_ * im_)
    check("cut_edge: V is exactly 1e4 and every |F| (= sqrt(Re^2 + Im^2)) equals 1e-4 exactly", vol == 1e4 and bool((absF == CUT).all()),
          f"{vol} {absF.min()} {absF.max()}")
    out["cut_edge"] = {"cell": list(ce_cell), "volume": vol, "site": {"frac": [0, 0, 0], "occupancy": 1.0},
                       "form_factor": "constant 1.0 for any Z and g2 (the kit's formFactor parameter)",
                       "reflections_within_k_max": int(len(Hc)), "absF_each": CUT,
                       "note": "with this form factor every reflection has |F| == 1e-4 exactly, so '>=' for '>' keeps all of them and '>' keeps none"}
    print(f"\ncut_edge: {len(Hc)} reflections within 1.2, every |F| == 1e-4 exactly")

    # ---- X-ray ----
    print("\n--- pymatgen X-ray pins (lambda 0.5 A) ---")
    xr = []
    for n in ("Au", "Si", "Mg", "beta_double_prime", "Cu3Au-occ", "WS2"):
        xr.append(xray_block(n, SS[n], SS[n]["pmg"]))
    out["xray"] = xr
    out["xray_worst_relative_difference"] = {x["name"]: x["worst_relative_difference"] for x in xr}

    # ---- structures ----
    out["structure_order"] = order
    out["structures"] = [SS[n]["entry"] for n in order]
    out["predictions"] = [{"name": n, "predicted": p, "measured": m_, "equal": p == m_} for n, p, m_ in PREDICTIONS]
    out["truth"] = truth

    # ---- write ----
    os.makedirs(OUT_DIR, exist_ok=True)
    dump = lambda o: json.dumps(o, separators=(",", ":"), allow_nan=False, default=float)
    with open(OUT, "w", encoding="utf-8") as f:
        keys = list(out.keys())
        f.write("{")
        for i, k in enumerate(keys):
            if k == "structures":
                f.write(f'"structures":[\n' + ",\n".join(dump(s) for s in out["structures"]) + "]")
            elif k == "xray":
                f.write('"xray":[\n' + ",\n".join(dump(s) for s in out["xray"]) + "]")
            elif k == "table":
                f.write('"table":[\n' + ",\n".join(dump(r) for r in out["table"]) + "]")
            elif k == "fe_samples":
                f.write('"fe_samples":{"g2":' + dump(out["fe_samples"]["g2"]) + ',"note":' + dump(out["fe_samples"]["note"])
                        + ',"rows":[\n' + ",\n".join(dump(r) for r in out["fe_samples"]["rows"]) + "]}")
            else:
                f.write(f'"{k}":' + dump(out[k]))
            f.write(",\n" if i < len(keys) - 1 else "\n")
        f.write("}\n")
    size = os.path.getsize(OUT)
    json.load(open(OUT, encoding="utf-8"))          # it parses
    print(f"\nwrote {OUT}  ({size / 1e6:.2f} MB, {len(out['structures'])} structures, {len(out['xray'])} X-ray sets)")
    print(f"{len(CHECKS)} self-checks ran, all passed (a failure would have stopped the run)")
    bad = [(n, p, m_) for n, p, m_ in PREDICTIONS if p != m_]
    print(f"{len(PREDICTIONS)} registered predictions, {len(PREDICTIONS) - len(bad)} equal, {len(bad)} differ")
    for n, p, m_ in bad:
        print(f"!!! PREDICTION DIFFERS: {n}: registered {p}, measured {m_}")
    if bad:
        sys.exit(3)


if __name__ == "__main__":
    main()
