"""Pins for the Hall-symbol reader and the Hall-only CIF.

Run: ../References/python/pymatgen/bin/python -I -B tools/pins/cif-symmetry/hall-pins.py
Versions: spglib 2.8.0, pymatgen-core 2026.10.16 (tools/pins/requirements/pymatgen.txt).

Writes Tests/DiffractoCoreTests/Fixtures/hall-symbol-pins.json:
  groups: for each of the 530 Hall numbers spglib knows, its Hall symbol, its IT number and every operation of the group in
          that setting (spglib.get_symmetry_from_database, centring translations included). One operation is 12
          characters: the rotation row by row as "-", "0", "+" and the translation in twelfths as three hex digits, each
          in 0 ... 11. The kit must build exactly this set from the symbol text.
  cifs:   planted CIFs that name their group by Hall symbol alone, with the atom count spglib's operations give (the
          truth) and the atom count pymatgen's CifParser builds (it does not read Hall symbols: with no operation loop it
          falls back to P1 with a warning; with a partial loop it expands with the loop as given, silently).
"""
import json, pathlib, warnings
import numpy as np
import spglib
from pymatgen.io.cif import CifParser

out = pathlib.Path(__file__).resolve().parents[3] / "Tests/DiffractoCoreTests/Fixtures/hall-symbol-pins.json"


def encode(r, t):
    rot = "".join("-0+"[int(v) + 1] for v in np.asarray(r).ravel())
    scaled = np.asarray(t, float) % 1.0 * 12
    assert np.allclose(scaled, np.rint(scaled), atol=1e-9), t
    tw = np.rint(scaled).astype(int) % 12
    return rot + "".join("0123456789ab"[int(v)] for v in tw)


groups = []
for h in range(1, 531):
    t = spglib.get_spacegroup_type(h)
    d = spglib.get_symmetry_from_database(h)
    ops = sorted({encode(r, tr) for r, tr in zip(d["rotations"], d["translations"])})
    assert len(ops) == len(d["rotations"]), h
    groups.append({"hall": h, "number": t.number, "symbol": t.hall_symbol, "operations": ops})


def hall_number(symbol):
    return next(g["hall"] for g in groups if g["symbol"] == symbol)


def xyz(r, t):
    parts = []
    for row, tt in zip(r, t):
        s = ""
        for c, v in zip("xyz", row):
            if v == 1: s += ("+" if s else "") + c
            elif v == -1: s += "-" + c
        tw = int(round(float(tt) % 1.0 * 12)) % 12
        if tw:
            from fractions import Fraction
            f = Fraction(tw, 12)
            s += f"+{f.numerator}/{f.denominator}"
        parts.append(s)
    return ",".join(parts)


def truth(symbol, sites):
    d = spglib.get_symmetry_from_database(hall_number(symbol))
    pts = []
    for el, f in sites:
        f = np.array(f, float)
        for r, t in zip(d["rotations"], d["translations"]):
            p = (r @ f + t) % 1.0
            if not any(e == el and np.allclose(((q - p + 0.5) % 1.0) - 0.5, 0, atol=1e-4) for e, q in pts):
                pts.append((el, p))
    return len(pts)


def cif(name, tag, symbol, cell, ops, sites):
    lines = [f"data_{name}", f"{tag} '{symbol}'"]
    for k, v in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"], cell):
        lines.append(f"_cell_{k} {v}")
    if ops is not None:
        lines += ["loop_", "_symmetry_equiv_pos_as_xyz"] + [f"'{o}'" for o in ops]
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
              "_atom_site_fract_z"]
    for i, (el, f) in enumerate(sites):
        lines.append(f"{el}{i + 1} {el} {f[0]} {f[1]} {f[2]}")
    return "\n".join(lines) + "\n"


def pymatgen_atoms(text):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = CifParser.from_str(text)
        return len(p.parse_structures(primitive=False, symmetrized=False)[0]), list(p.warnings)


cubic_point = [xyz(r, (0, 0, 0)) for r in spglib.get_symmetry_from_database(hall_number("-P 4 2 3"))["rotations"]]
nacl = [("Na", (0, 0, 0)), ("Cl", (0.5, 0.5, 0.5))]
specs = [
    ("rock salt, Hall only, no loop", "_space_group_name_Hall", "-F 4 2 3", (5.64, 5.64, 5.64, 90, 90, 90), None, nacl),
    ("rock salt, legacy Hall tag, no loop", "_symmetry_space_group_name_Hall", "-F 4 2 3",
     (5.64, 5.64, 5.64, 90, 90, 90), None, nacl),
    ("P 21 general site, Hall only, no loop", "_space_group_name_Hall", "P 2yb", (5, 6, 7, 90, 100, 90), None,
     [("Si", (0.12, 0.23, 0.34))]),
    ("P 31 2 1 with origin shift, Hall only, no loop", "_space_group_name_Hall", "P 31 2 (0 0 4)",
     (4.9, 4.9, 5.4, 90, 90, 120), None, [("Si", (0.47, 0.0, 0.0)), ("O", (0.41, 0.27, 0.12))]),
    ("rhombohedral axes, Hall only, no loop", "_space_group_name_Hall", "-P 3* 2", (5.1, 5.1, 5.1, 55, 55, 55), None,
     [("Al", (0.35, 0.35, 0.35)), ("O", (0.55, 0.95, 0.25))]),
    ("rock salt, Hall only, 48 point operations", "_space_group_name_Hall", "-F 4 2 3", (5.64, 5.64, 5.64, 90, 90, 90),
     cubic_point, nacl),
]
cifs = []
for name, tag, symbol, cell, ops, sites in specs:
    text = cif("t", tag, symbol, cell, ops, sites)
    n, w = pymatgen_atoms(text)
    order = len(next(g for g in groups if g["symbol"] == symbol)["operations"])
    cifs.append({"name": name, "text": text, "hall": symbol, "truthAtoms": truth(symbol, sites), "pymatgenAtoms": n,
                 "pymatgenWarnings": w, "listed": None if ops is None else len(ops), "groupOrder": order})

json.dump({"spglib": spglib.__version__, "pymatgen-core": "2026.10.16", "groups": groups, "cifs": cifs},
          open(out, "w"), separators=(",", ":"))
print(out, len(groups), sum(len(g["operations"]) for g in groups))
for c in cifs:
    print(c["name"], c["truthAtoms"], c["pymatgenAtoms"], c["pymatgenWarnings"])
