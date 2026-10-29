"""DiffractoKit 0.0.4: the tetragonal pins, in three axis settings.

Run, from any directory, with the pinned oracle:

    ../References/python/orix/bin/python -I -B tools/kit-0.0.4/tetragonal-pins.py

It writes Tests/DiffractoCoreTests/Fixtures/tetragonal-pins.json BEFORE the Swift fix exists. Oracle: pymatgen 2026.9.24
(pymatgen-core 2026.9.23) SpacegroupAnalyzer on spglib 2.8.0, NumPy 2.5.4.

Tolerance: symprec = 1e-4 (Angstrom, spglib's), as slice 1's pins (tools/slice-1/pins.py SYMPREC). The inputs are
exact (coordinates written to 4 decimals and expanded by the space group's own operations, so every equivalent site is
exact to 1e-15): the answer does not depend on symprec between 1e-6 and 1e-2, which the script checks (symprec 1e-2 gives the
same space group number). The kit's own position tolerance is a fixed 5e-3 fractional (SymmetryDetection.positionTolerance);
the metric tolerance is 1e-9 relative.

Structures (published, conventional cells):
  rutile TiO2   P4_2/mnm (136)  a = 4.5937, c = 2.9587 A; Ti 2a (0,0,0), O 4f (u,u,0) with u = 0.3048
                (Burdett, Hughbanks, Miller, Richardson, Smith, J. Am. Chem. Soc. 109, 3639 (1987), neutron data, 15 K).
  anatase TiO2  I4_1/amd (141, origin choice 1)  a = 3.7845, c = 9.5143 A; Ti 4a (0,0,0), O 8e (0,0,z) with z = 0.2081
                (Horn, Schwebdtfeger, Meyer, Z. Kristallogr. 136, 273 (1972); Ti-O 1.934 A is checked below).
  beta-Sn       I4_1/amd (141, origin choice 1)  a = 5.8318, c = 3.1818 A; Sn 4a (0,0,0)
                (Wyckoff, Crystal Structures; Sn-Sn 3.022 A is checked below).
Controls:
  ortho_a_equals_c  an orthorhombic cell with a = c by accident (4.0, 5.0, 4.0 A) holding Ti at the origin and O on the
                    general position (0.3, 0.2, 0.1) of Pmmm (eight O): mmm, not 4-fold in any setting.
  cubic             SrTiO3 Pm-3m, a = 3.905 A (Sr (0,0,0), Ti (1/2,1/2,1/2), O (1/2,1/2,0)): m-3m.

Settings: "c" has the 4-fold along c (the conventional cell); "a" has it along a (new a, b, c = old c, a, b: lengths
(c, a, a)); "b" has it along b (new a, b, c = old b, c, a: lengths (a, c, a)). Both are cyclic permutations, so the cell
stays right-handed; a site (x, y, z) of the "c" setting is (z, x, y) in "a" and (y, z, x) in "b". The
controls take the same three permutations: the orthorhombic control is a = c in its "c" listing, a = b in "a" and b = c
in "b" (so each setting puts its accidental equal pair on a different pair of axes, and the kit must not read 4-fold in
any), and the cubic one is invariant.

Per entry the fixture carries: the cell, the sites, pymatgen's space-group number, point-group symbol, the Laue class
(the 32 -> 11 table of slice 1), the number of point operations of the point group (get_point_group_operations), the number of
space-group operations on the cell as given (get_symmetry_operations; centring translations counted, lattice
translations not), and the order of the Laue group (2 x the point-group's proper part: the count of {+-R}).
"""
# Run under the suite's memory cap (a process cannot cap itself on macOS 27):
#   ../References/python/orix/bin/python -I -B tools/kit-0.0.4/tetragonal-pins.py
import json
import math
import os
import warnings
from importlib.metadata import version as dist_version

import numpy as np
import spglib

import pymatgen.core
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(KIT, "Tests", "DiffractoCoreTests", "Fixtures", "tetragonal-pins.json")
SYMPREC = 1e-4

LAUE_OF_POINT_GROUP = {
    "1": "-1", "-1": "-1", "2": "2/m", "m": "2/m", "2/m": "2/m", "222": "mmm", "mm2": "mmm", "mmm": "mmm",
    "4": "4/m", "-4": "4/m", "4/m": "4/m", "422": "4/mmm", "4mm": "4/mmm", "-42m": "4/mmm", "-4m2": "4/mmm",
    "4/mmm": "4/mmm", "3": "-3", "-3": "-3", "32": "-3m", "3m": "-3m", "-3m": "-3m",
    "6": "6/m", "-6": "6/m", "6/m": "6/m", "622": "6/mmm", "6mm": "6/mmm", "-6m2": "6/mmm", "-62m": "6/mmm",
    "6/mmm": "6/mmm", "23": "m-3", "m-3": "m-3", "432": "m-3m", "-43m": "m-3m", "m-3m": "m-3m",
}
LAUE_ORDER = {"-1": 2, "2/m": 4, "mmm": 8, "4/m": 8, "4/mmm": 16, "-3": 6, "-3m": 12, "6/m": 12, "6/mmm": 24,
              "m-3": 24, "m-3m": 48}


def wrap(x):
    return [float(v % 1.0) if abs(v % 1.0 - 1.0) > 1e-12 else 0.0 for v in x]


def permute(setting, lengths, coords):
    """lengths (a, b, c) and fractional coordinates of the c-setting, re-expressed in setting 'c', 'a' or 'b'."""
    a, b, c = lengths
    if setting == "c":
        return (a, b, c), [list(x) for x in coords]
    if setting == "a":      # new axes (a, b, c) = old (c, a, b): site (x, y, z) -> (z, x, y)
        return (c, a, b), [[x[2], x[0], x[1]] for x in coords]
    if setting == "b":      # new axes = old (b, c, a): site (x, y, z) -> (y, z, x)
        return (b, c, a), [[x[1], x[2], x[0]] for x in coords]
    raise ValueError(setting)


def build(sg, lengths, species, coords):
    return Structure.from_spacegroup(sg, Lattice.orthorhombic(*lengths), species, coords)


def sites_of(structure):
    return [{"z": int(site.specie.Z), "frac": wrap(site.frac_coords)} for site in structure]


def pin_of(name, setting, lengths, species_z, frac):
    """pymatgen's analysis of the cell given by lengths (90 deg angles) and the sites (z, frac)."""
    lattice = Lattice.orthorhombic(*lengths)
    structure = Structure(lattice, [int(z) for z in species_z], frac, coords_are_cartesian=False)
    sga = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    number = sga.get_space_group_number()
    # the answer must not hang on the symprec: the inputs are exact
    for other in (1e-2, 1e-3):
        assert SpacegroupAnalyzer(structure, symprec=other).get_space_group_number() == number, (name, setting, other)
    pg = sga.get_point_group_symbol()
    laue = LAUE_OF_POINT_GROUP[pg]
    point_ops = len(sga.get_point_group_operations())
    space_ops = len(sga.get_symmetry_operations())
    # point operations of the point group x inversion if absent = the Laue group's order
    assert LAUE_ORDER[laue] in (point_ops, 2 * point_ops), (name, point_ops, laue)
    return {
        "name": name, "setting": setting,
        "cell": {"a": lengths[0], "b": lengths[1], "c": lengths[2], "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "sites": [{"z": int(z), "frac": wrap(f)} for z, f in zip(species_z, frac)],
        "space_group_number": int(number), "point_group": pg, "laue_class": laue,
        "point_group_operations": int(point_ops), "laue_group_order": LAUE_ORDER[laue],
        "space_group_operations": int(space_ops),
    }


def main():
    warnings.filterwarnings("ignore")
    entries = []
    base = {}

    rutile = build("P4_2/mnm", (4.5937, 4.5937, 2.9587), ["Ti", "O"], [[0, 0, 0], [0.3048, 0.3048, 0]])
    anatase = build("I4_1/amd", (3.7845, 3.7845, 9.5143), ["Ti", "O"], [[0, 0, 0], [0, 0, 0.2081]])
    beta_sn = build("I4_1/amd", (5.8318, 5.8318, 3.1818), ["Sn"], [[0, 0, 0]])
    assert len(rutile) == 6 and len(anatase) == 12 and len(beta_sn) == 4
    ti_o = min(anatase.get_distance(0, j) for j in range(len(anatase)) if anatase[j].specie.Z == 8)
    assert abs(ti_o - 1.934) < 0.01, ti_o
    sn_sn = min(beta_sn.get_distance(0, j) for j in range(1, len(beta_sn)))
    assert abs(sn_sn - 3.022) < 0.005, sn_sn
    ti_o_rutile = min(rutile.get_distance(0, j) for j in range(len(rutile)) if rutile[j].specie.Z == 8)
    assert abs(ti_o_rutile - 1.946) < 0.01, ti_o_rutile      # rutile's apical Ti-O 1.946 A

    ortho = build("Pmmm", (4.0, 5.0, 4.0), ["Ti", "O"], [[0, 0, 0], [0.3, 0.2, 0.1]])
    cubic = build("Pm-3m", (3.905, 3.905, 3.905), ["Sr", "Ti", "O"], [[0, 0, 0], [.5, .5, .5], [.5, .5, 0]])
    assert len(ortho) == 9 and len(cubic) == 5

    sources = [
        ("rutile_TiO2", rutile, (4.5937, 4.5937, 2.9587)),
        ("anatase_TiO2", anatase, (3.7845, 3.7845, 9.5143)),
        ("beta_Sn", beta_sn, (5.8318, 5.8318, 3.1818)),
        ("ortho_a_equals_c", ortho, (4.0, 5.0, 4.0)),
        ("cubic_SrTiO3", cubic, (3.905, 3.905, 3.905)),
    ]
    for name, structure, lengths in sources:
        zs = [int(s.specie.Z) for s in structure]
        coords = [list(s.frac_coords) for s in structure]
        for setting in ("c", "a", "b"):
            L, F = permute(setting, lengths, coords)
            entries.append(pin_of(name, setting, L, zs, F))

    # self-checks: the tetragonal structures and the cubic one have the same group in every setting, the control is mmm
    by = {}
    for e in entries:
        by.setdefault(e["name"], []).append(e)
    for name, group in by.items():
        keys = {(e["space_group_number"], e["point_group"], e["laue_class"], e["point_group_operations"],
                 e["space_group_operations"]) for e in group}
        assert len(keys) == 1, (name, keys)
    assert {e["laue_class"] for e in by["rutile_TiO2"] + by["anatase_TiO2"] + by["beta_Sn"]} == {"4/mmm"}
    assert {e["laue_class"] for e in by["ortho_a_equals_c"]} == {"mmm"}
    assert {e["laue_class"] for e in by["cubic_SrTiO3"]} == {"m-3m"}

    out = {
        "versions": {"pymatgen_core": dist_version("pymatgen-core"), "pymatgen": dist_version("pymatgen"),
                     "spglib": spglib.__version__, "numpy": np.__version__,
                     "python": __import__("sys").version.split()[0]},
        "symprec": SYMPREC, "k_max": 1.2,
        "entries": entries,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    for e in entries:
        print("%-18s %s  sg %3d  pg %-6s laue %-5s  point-ops %2d  space-ops %3d  sites %2d"
              % (e["name"], e["setting"], e["space_group_number"], e["point_group"], e["laue_class"],
                 e["point_group_operations"], e["space_group_operations"], len(e["sites"])))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
