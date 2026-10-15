"""DiffractoKit slice 1b: the pin generator (registration docs/registration-slice-1b.md, section 2).

Run, from any directory, with the pinned oracle, UNDER THE SUITE MEMORY WATCHDOG (owner, 2026-10-15 12:30: setrlimit and
`ulimit -v` fail on macOS 27, so the cap is the watchdog's, 8 GB, exit 137; one such worker at a time):

    <suite>/tools/memcap.sh -- References/pymatgen-venv/bin/python -I -B tools/slice-1b/pins.py

It writes Tests/DiffractoCoreTests/Fixtures/slice-1b-pins.json BEFORE any Swift for the reduction exists.

Oracles (all from the pinned venv): spglib 2.8.0 `niggli_reduce` and `get_symmetry_dataset` / `standardize_cell`;
pymatgen-core 2026.9.23 `Lattice.get_niggli_reduced_lattice` as the second, independent Niggli implementation. The ring
list in the conventional setting uses slice 1's own definitions (tools/slice-1/pins.py, imported, not copied) on spglib's
standardised cell. Agreement is agreement, not truth: truth (the Niggli conditions, P unimodular, the group-theory
counts) is checked in the Swift tests directly.

Inputs: slice-1-pins.json (cell and expanded sites of every slice-1 structure), the named non-conventional cells built
below, and RANDOM_PER_LATTICE random cells of each of the 14 Bravais lattices under a random unimodular transformation
(|entries| <= 2), seed SEED.

Conventions of the pins (the Swift tests read them): a crystal is its six cell parameters and its fractional sites; the
lattice matrix every oracle is given is pymatgen's `Lattice.from_parameters(..., vesta=True)` (a along x, b in xy), the
same lattice `UnitCell` builds. `bravais` is the lattice type of spglib's LATTICE-ONLY dataset (one atom at the origin):
the holohedry's space group (Pm-3m cP, Im-3m cI, Fm-3m cF, ...), whose centring letter and family give the name.
"""
import datetime
import importlib.util
import json
import math
import os
import sys
import warnings
from importlib.metadata import version as dist_version

import numpy as np
import spglib

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURES = os.path.join(KIT, "Tests", "DiffractoCoreTests", "Fixtures")
SLICE1_PINS = os.path.join(FIXTURES, "slice-1-pins.json")
OUT = os.path.join(FIXTURES, "slice-1b-pins.json")

SEED = 20261015
RANDOM_PER_LATTICE = 200      # registration 2: 200; the report states it if this is reduced
SYMPREC = 1e-4
NIGGLI_TOL = 1e-5             # spglib niggli_reduce(eps=) and pymatgen get_niggli_reduced_lattice(tol=)
K_MAX = 1.2
NIGGLI_AGREE = 1e-8           # the two oracles must agree on the reduced parameters (relative lengths, degrees)

# slice 1's generator, imported as a module (its main() is guarded)
_spec = importlib.util.spec_from_file_location("slice1_pins", os.path.join(HERE, "..", "slice-1", "pins.py"))
s1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s1)

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    assert ok, f"SELF-CHECK FAILED: {name} {detail}"


# ----------------------------------------------------------------------------------------------------------------
# Lattice type of a space-group number (the lattice-only dataset is the holohedry's, so this names the LATTICE)
# ----------------------------------------------------------------------------------------------------------------
def bravais_of(number, international):
    centring = international.strip()[0]
    if number <= 2:
        fam = "a"
    elif number <= 15:
        fam = "m"
    elif number <= 74:
        fam = "o"
    elif number <= 142:
        fam = "t"
    elif number <= 194:
        fam = "h"
    else:
        fam = "c"
    if fam == "h" and centring == "R":
        return "hR"
    return fam + centring


BRAVAIS_NAMES = ["aP", "mP", "mC", "oP", "oC", "oI", "oF", "tP", "tI", "hR", "hP", "cP", "cI", "cF"]
HOLOHEDRY_ORDER = {"aP": 2, "mP": 4, "mC": 4, "oP": 8, "oC": 8, "oI": 8, "oF": 8, "tP": 16, "tI": 16, "hR": 12,
                   "hP": 24, "cP": 48, "cI": 48, "cF": 48}
HEX_AXES_CLASSES = ("-3", "-3m", "6/m", "6/mmm")


def params_of(matrix):
    L = Lattice(matrix)
    a, b, c = L.abc
    al, be, ga = L.angles
    return [float(a), float(b), float(c), float(al), float(be), float(ga)]


def lattice_matrix(params):
    return np.array(Lattice.from_parameters(*params, vesta=True).matrix)


def reduced_params_spglib(A):
    return params_of(np.array(spglib.niggli_reduce(np.array(A), eps=NIGGLI_TOL)))


def reduced_params_pymatgen(A):
    return [float(x) for x in Lattice(A).get_niggli_reduced_lattice(tol=NIGGLI_TOL).parameters]


def params_close(p, q, rel=NIGGLI_AGREE, deg=NIGGLI_AGREE):
    return all(abs(p[i] - q[i]) <= rel * max(abs(p[i]), 1.0) for i in range(3)) and \
        all(abs(p[i] - q[i]) <= deg for i in range(3, 6))


def worst_param_diff(p, q):
    return max(max(abs(p[i] - q[i]) / max(abs(p[i]), 1.0) for i in range(3)), max(abs(p[i] - q[i]) for i in range(3, 6)))


def lattice_only(A):
    """spglib's dataset of one atom at the origin: the lattice's holohedry (Bravais name, order, conventional cell)."""
    ds = spglib.get_symmetry_dataset((np.array(A), [[0.0, 0.0, 0.0]], [1]), symprec=SYMPREC)
    std = np.array(ds.std_lattice)
    lens = sorted(float(x) for x in np.linalg.norm(std, axis=1))
    return {
        "bravais": bravais_of(int(ds.number), ds.international),
        "holohedry_space_group": int(ds.number),
        "holohedry_order": int(len(ds.rotations)),
        "conv_volume": float(abs(np.linalg.det(std))),
        "conv_params": params_of(std),
        "conv_lengths_sorted": lens,
        "conv_unique_b": float(np.linalg.norm(std[1])),
    }


# ----------------------------------------------------------------------------------------------------------------
# Structures: oracle entry for (params, Z, occ, frac)
# ----------------------------------------------------------------------------------------------------------------
def structure_entry(name, source, params, Z, occ, frac):
    A = lattice_matrix(params)
    V = abs(float(np.linalg.det(A)))
    sp = reduced_params_spglib(A)
    pm = reduced_params_pymatgen(A)
    check(f"{name}: spglib and pymatgen reduce to the same parameters (worst {worst_param_diff(sp, pm):.2e})",
          params_close(sp, pm), f"{sp} {pm}")
    check(f"{name}: reduced cell keeps the volume", abs(abs(np.linalg.det(np.array(spglib.niggli_reduce(A, eps=NIGGLI_TOL)))) - V) <= 1e-9 * V)
    lat = lattice_only(A)
    frac = np.array(frac, dtype=float)
    Zs = [int(z) for z in Z]
    ds = spglib.get_symmetry_dataset((A, frac, Zs), symprec=SYMPREC)
    pg = ds.pointgroup
    check(f"{name}: point group {pg!r} in the 32 -> 11 table", pg in s1.LAUE_OF_POINT_GROUP)
    laue_crystal = s1.LAUE_OF_POINT_GROUP[pg]
    # the Laue class of the operations IN THIS CELL: a cell that is a supercell cannot hold all of the crystal's
    # rotations (a 1x1x7 supercell of a cubic crystal holds 4/mmm), and the kit's detection works inside the cell
    rots = np.unique(np.array(ds.rotations).reshape(-1, 3, 3), axis=0)
    laue = s1.laue_class_from_group(s1.laue_rotations(rots))
    n_rot = int(len(ds.rotations))
    ratio = lat["conv_volume"] / V
    check(f"{name}: conventional-cell volume ratio {ratio} is an integer", abs(ratio - round(ratio)) < 1e-6)
    conv_ops = n_rot * int(round(ratio))
    # the standardised (conventional, not primitive) cell and its rings: slice 1's own definitions on it
    std = spglib.standardize_cell((A, frac, Zs), to_primitive=False, no_idealize=False, symprec=SYMPREC)
    check(f"{name}: spglib standardises the cell", std is not None)
    sL, sP, sZ = std
    sparams = params_of(np.array(sL))
    st = Structure(Lattice(np.array(sL)), [int(z) for z in sZ], np.array(sP), coords_are_cartesian=False)
    sga, Rs, ts, _ = s1.operations_of(st)
    LR = s1.laue_rotations(Rs)
    laue_std = s1.laue_class_from_group(LR)
    check(f"{name}: Laue class of the standardised cell's group ({laue_std}) equals the table's ({laue_crystal})",
          laue_std == laue_crystal)
    Ls, As, Gs, Vs, Bs = s1.lattice_block(*sparams)
    metric_hex = abs(sparams[0] - sparams[1]) < 1e-9 * sparams[0] and abs(sparams[5] - 120) < 1e-6 \
        and abs(sparams[3] - 90) < 1e-6 and abs(sparams[4] - 90) < 1e-6
    hex_axes = (laue_std in HEX_AXES_CLASSES) and metric_hex
    cache = s1.OrbitCache(LR)
    r = s1.reflections_for(Ls, As, Bs, K_MAX, Rs, ts, LR, hex_axes, cache, detail=False)
    fams = [{"g": f["g"], "multiplicity": f["multiplicity"], "name": f["name"], "orbits": f["orbits"]} for f in r["families"]]
    print(f"  {name:<34s} {lat['bravais']}  laue {laue:<6s}({laue_crystal:<6s}) ops(file) {n_rot:4d} ops(conv) {conv_ops:4d}  "
          f"reduced {['%.4f' % x for x in sp]}  rings {len(fams)}")
    return {
        "name": name, "source": source,
        "cell": dict(zip(["a", "b", "c", "alpha", "beta", "gamma"], params)),
        "sites": [{"z": Zs[i], "frac": [float(x) for x in frac[i]], "occupancy": float(occ[i])} for i in range(len(Zs))],
        "niggli_spglib": sp, "niggli_pymatgen": pm,
        "lattice": lat,
        "laue_class": laue, "crystal_laue_class": laue_crystal, "point_group": pg, "space_group_number": int(ds.number),
        "file_operations_count": n_rot, "conventional_operations_count": conv_ops,
        "conv_std_params": sparams,
        "rings": {"k_max": K_MAX, "count": r["count"], "allowed_count": r["allowed_count"], "families": fams},
    }


# ----------------------------------------------------------------------------------------------------------------
# The named non-conventional cells (registration 2)
# ----------------------------------------------------------------------------------------------------------------
def transformed(A, frac, P):
    """P: columns are the new basis vectors in terms of the old (x_old = P x_new). Returns params, frac_new."""
    P = np.array(P, dtype=float)
    A2 = P.T @ A
    Pinv = np.linalg.inv(P)
    f2 = (Pinv @ np.array(frac).T).T
    f2 = f2 - np.floor(f2 + 1e-9)
    return params_of(A2), f2


def named_cells(slice1):
    out = []
    a = 4.0782
    A = a / 2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    out.append(("au_fcc_primitive", params_of(A), [79], [1.0], [[0, 0, 0]]))
    a = 2.8665
    A = a / 2 * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=float)
    out.append(("fe_bcc_primitive", params_of(A), [26], [1.0], [[0, 0, 0]]))
    # beta'' C-centred primitive cell: spglib's primitive of the 22-atom cell
    b = next(s for s in slice1["structures"] if s["name"] == "beta_double_prime")
    cellb = [b["cell"][k] for k in ("a", "b", "c", "alpha", "beta", "gamma")]
    stb = Structure(Lattice.from_parameters(*cellb, vesta=True), [s["z"] for s in b["sites"]],
                    [s["frac"] for s in b["sites"]])
    prim = SpacegroupAnalyzer(stb, symprec=SYMPREC).find_primitive()
    check("beta'' primitive cell has 11 atoms", len(prim) == 11, f"{len(prim)}")
    out.append(("beta_double_prime_primitive", params_of(prim.lattice.matrix), [int(s.specie.Z) for s in prim],
                [1.0] * len(prim), [list(map(float, s.frac_coords)) for s in prim]))
    # rutile with the 4-fold axis along a: (a', b', c') = (c, a, b)
    ar, cr, u = 4.5937, 2.9587, 0.3053
    sites = [(22, (0, 0, 0)), (22, (.5, .5, .5)), (8, (u, u, 0)), (8, (-u, -u, 0)), (8, (.5 + u, .5 - u, .5)),
             (8, (.5 - u, .5 + u, .5))]
    out.append(("rutile_unique_axis_a", [cr, ar, ar, 90.0, 90.0, 90.0], [s[0] for s in sites], [1.0] * 6,
                [[z, x, y] for _, (x, y, z) in sites]))
    # magnesium with gamma = 60 degrees: a2' = a1 + a2
    am, cm = 3.2094, 5.2108
    Ah = lattice_matrix([am, am, cm, 90, 90, 120])
    fh = [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    p, f = transformed(Ah, fh, [[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    out.append(("magnesium_gamma_60", p, [12, 12], [1.0, 1.0], f.tolist()))
    # corundum in rhombohedral axes: primitive cell of R-3c (hexagonal axes a=4.7602, c=12.9933)
    cor = Structure.from_spacegroup("R-3c", Lattice.from_parameters(4.7602, 4.7602, 12.9933, 90, 90, 120, vesta=True),
                                    ["Al", "O"], [[0, 0, 0.3523], [0.3064, 0, 0.25]])
    check("corundum hexagonal cell holds 30 atoms", len(cor) == 30, f"{len(cor)}")
    prim = SpacegroupAnalyzer(cor, symprec=SYMPREC).find_primitive()
    check("corundum rhombohedral cell holds 10 atoms", len(prim) == 10, f"{len(prim)}")
    out.append(("corundum_rhombohedral_axes", params_of(prim.lattice.matrix), [int(s.specie.Z) for s in prim],
                [1.0] * len(prim), [list(map(float, s.frac_coords)) for s in prim]))
    # the 7x supercell fixture of Gate B (a P1 list)
    path = os.path.join(FIXTURES, "refuter_F3_super7.cif")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p7 = CifParser(path, frac_tolerance=0.0)
        st7 = p7.parse_structures(primitive=False)[0]
    Z, occ, frac = s1.sites_of(st7)
    out.append(("supercell_7x", [float(x) for x in st7.lattice.parameters], [int(z) for z in Z], [float(o) for o in occ],
                [list(map(float, f)) for f in frac]))
    # a triclinic cell whose Niggli-reduced form is hexagonal: hcp Mg in a skewed unimodular basis
    M = np.array([[1, 2, 0], [0, 1, 0], [1, 1, 1]])
    check("skew matrix is unimodular", abs(round(np.linalg.det(M))) == 1)
    p, f = transformed(Ah, fh, M)
    out.append(("triclinic_hexagonal_reduced", p, [12, 12], [1.0, 1.0], f.tolist()))
    return out


# ----------------------------------------------------------------------------------------------------------------
# Random cells of every Bravais lattice under random unimodular transformations
# ----------------------------------------------------------------------------------------------------------------
PRIMITIVE_BASIS = {   # rows: primitive vectors in units of the conventional a, b, c
    "P": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "C": [[.5, .5, 0], [-.5, .5, 0], [0, 0, 1]],
    "I": [[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]],
    "F": [[0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
    "R": [[2 / 3, 1 / 3, 1 / 3], [-1 / 3, 1 / 3, 1 / 3], [-1 / 3, -2 / 3, 1 / 3]],
}


def conventional_params(kind, rng):
    u = rng.uniform
    if kind == "aP":
        return [u(3, 9), u(3, 9), u(3, 9), u(60, 118), u(60, 118), u(60, 118)]
    if kind in ("mP", "mC"):
        return [u(3, 9), u(3, 9), u(3, 9), 90.0, u(95, 125), 90.0]
    if kind in ("oP", "oC", "oI", "oF"):
        return [u(3, 9), u(3, 9), u(3, 9), 90.0, 90.0, 90.0]
    if kind in ("tP", "tI"):
        a = u(3, 8)
        return [a, a, a * u(0.5, 2.2), 90.0, 90.0, 90.0]
    if kind == "hP":
        a = u(3, 7)
        return [a, a, a * u(0.8, 3.0), 90.0, 90.0, 120.0]
    if kind == "hR":
        a = u(3, 6)
        return [a, a, a * u(2.0, 4.5), 90.0, 90.0, 120.0]
    a = u(3, 8)
    return [a, a, a, 90.0, 90.0, 90.0]


def random_unimodular(rng):
    while True:
        M = rng.integers(-2, 3, size=(3, 3))
        if abs(round(np.linalg.det(M))) == 1:
            return M


def random_cells(rng):
    cells = []
    redraws = 0
    worst_oracle = 0.0
    for kind in BRAVAIS_NAMES:
        got = 0
        while got < RANDOM_PER_LATTICE:
            conv = conventional_params(kind, rng)
            Ac = lattice_matrix(conv)
            Ap = np.array(PRIMITIVE_BASIS[kind[1]]) @ Ac
            M = random_unimodular(rng)
            Ain = M @ Ap
            if np.linalg.norm(Ain, axis=1).max() > 60:
                redraws += 1
                continue
            lat = lattice_only(Ain)
            if lat["bravais"] != kind:
                redraws += 1            # a special ratio made the lattice something else: the draw is not this lattice
                continue
            sp = reduced_params_spglib(Ain)
            pm = reduced_params_pymatgen(Ain)
            if not params_close(sp, pm):
                redraws += 1            # the two oracles' epsilon decisions differ here: no agreed reference (counted)
                continue
            worst_oracle = max(worst_oracle, worst_param_diff(sp, pm))
            p_in = params_of(Ain)
            cells.append({"intended": kind, "cell": dict(zip(["a", "b", "c", "alpha", "beta", "gamma"], p_in)),
                          "niggli_spglib": sp, "niggli_pymatgen": pm, "lattice": lat, "unimodular": M.tolist()})
            got += 1
    return cells, redraws, worst_oracle


def main():
    versions = {
        "python": sys.version, "pymatgen": dist_version("pymatgen"), "pymatgen_core": dist_version("pymatgen-core"),
        "spglib": spglib.__version__, "numpy": np.__version__,
    }
    print("oracle:", json.dumps(versions))
    with open(SLICE1_PINS, encoding="utf-8") as f:
        slice1 = json.load(f)
    structures = []
    print("\n--- slice 1's fixtures (P1 sites, as the Swift tests build them) ---")
    for s in slice1["structures"]:
        sites = s.get("sites_as_written") or s["sites"]
        params = [s["cell"][k] for k in ("a", "b", "c", "alpha", "beta", "gamma")]
        e = structure_entry(s["name"], "slice-1", params, [x["z"] for x in sites], [x["occupancy"] for x in sites],
                            [x["frac"] for x in sites])
        # slice 1's own pin for the same crystal agrees where the setting is unchanged
        if abs(e["lattice"]["conv_volume"] / abs(np.linalg.det(lattice_matrix(params))) - 1) < 1e-9:
            check(f"{s['name']}: conventional operation count equals slice 1's pin",
                  e["conventional_operations_count"] == len(s["operations"]),
                  f"{e['conventional_operations_count']} vs {len(s['operations'])}")
        structures.append(e)
    print("\n--- named non-conventional cells ---")
    for name, params, Z, occ, frac in named_cells(slice1):
        structures.append(structure_entry(name, "named", params, Z, occ, frac))

    # registration 3: the counts the registration states
    by = {s["name"]: s for s in structures}
    check("truth: Au primitive is cF with 192 conventional operations, m-3m",
          by["au_fcc_primitive"]["lattice"]["bravais"] == "cF" and by["au_fcc_primitive"]["conventional_operations_count"] == 192
          and by["au_fcc_primitive"]["laue_class"] == "m-3m")
    check("truth: Fe primitive is cI with 96, m-3m", by["fe_bcc_primitive"]["lattice"]["bravais"] == "cI"
          and by["fe_bcc_primitive"]["conventional_operations_count"] == 96)
    check("truth: beta'' primitive is mC with 8, 2/m", by["beta_double_prime_primitive"]["lattice"]["bravais"] == "mC"
          and by["beta_double_prime_primitive"]["conventional_operations_count"] == 8
          and by["beta_double_prime_primitive"]["laue_class"] == "2/m")
    check("truth: rutile with a unique: 4/mmm, 16", by["rutile_unique_axis_a"]["laue_class"] == "4/mmm"
          and by["rutile_unique_axis_a"]["conventional_operations_count"] == 16)
    check("truth: Mg gamma 60: 6/mmm, 24", by["magnesium_gamma_60"]["laue_class"] == "6/mmm"
          and by["magnesium_gamma_60"]["conventional_operations_count"] == 24)
    check("truth: corundum rhombohedral axes: -3m, hR, 36 conventional", by["corundum_rhombohedral_axes"]["laue_class"] == "-3m"
          and by["corundum_rhombohedral_axes"]["lattice"]["bravais"] == "hR"
          and by["corundum_rhombohedral_axes"]["conventional_operations_count"] == 36)
    check("truth: triclinic-looking cell reduces to a hexagonal lattice",
          by["triclinic_hexagonal_reduced"]["lattice"]["bravais"] == "hP")

    rng = np.random.default_rng(SEED)
    print(f"\n--- random cells: {RANDOM_PER_LATTICE} per Bravais lattice, seed {SEED} ---")
    rcells, redraws, worst_oracle = random_cells(rng)
    print(f"{len(rcells)} cells, {redraws} draws rejected (a special ratio changed the lattice, or the two Niggli oracles "
          f"disagreed beyond {NIGGLI_AGREE}); worst spglib/pymatgen reduced-parameter difference accepted {worst_oracle:.2e}")

    out = {
        "versions": versions,
        "generated_at": datetime.datetime.now(datetime.UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "definitions": {
            "seed": SEED, "random_per_lattice": RANDOM_PER_LATTICE, "symprec": SYMPREC, "niggli_eps": NIGGLI_TOL,
            "k_max": K_MAX,
            "bravais": "name of spglib's lattice-only dataset (one atom at the origin): family letter + centring letter",
            "conventional_operations": "len(dataset.rotations) of the input cell times the integer ratio of the lattice's "
                                       "conventional-cell volume to the input cell's volume",
            "rings": "slice 1's definitions on spglib's standardised (conventional, idealised) cell, symmetry-allowed "
                     "reflections at k_max, names in that setting",
            "niggli_agreement": NIGGLI_AGREE,
        },
        "structures": structures, "random_cells": rcells,
        "random_summary": {"rejected_draws": redraws, "worst_oracle_difference": worst_oracle},
    }
    dump = lambda o: json.dumps(o, separators=(",", ":"), allow_nan=False)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write('{"versions":' + dump(out["versions"]) + ',\n')
        f.write('"generated_at":' + dump(out["generated_at"]) + ',\n')
        f.write('"definitions":' + dump(out["definitions"]) + ',\n')
        f.write('"random_summary":' + dump(out["random_summary"]) + ',\n')
        f.write('"structures":[\n' + ",\n".join(dump(s) for s in structures) + '],\n')
        f.write('"random_cells":[\n' + ",\n".join(dump(c) for c in rcells) + ']}\n')
    print(f"\nwrote {OUT} ({os.path.getsize(OUT) / 1e6:.2f} MB, {len(structures)} structures, {len(rcells)} random cells)")
    print(f"{len(CHECKS)} self-checks ran; failed: {[n for n, ok in CHECKS if not ok]}")


if __name__ == "__main__":
    main()
