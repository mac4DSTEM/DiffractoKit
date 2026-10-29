"""Pins for the CIF symmetry-completeness check (0.9.0).

Run: ../References/python/pymatgen/bin/python -I -B tools/pins/cif-symmetry/pins.py
Versions: spglib 2.8.0; the fixture's "pymatgen" key is a literal, 2026.9.24, the metapackage version. The code on disk in
References/python/pymatgen is pymatgen-core 2026.10.16 (checked 2026-10-21: regenerating gave a byte-identical file, so the numbers
are the same under both). Other fixtures record pymatgen-core 2026.9.23 from an earlier venv.

Writes Tests/DiffractoCoreTests/Fixtures/cif-symmetry-pins.json:
  symbols: every Hermann-Mauguin spelling spglib lists (short, full, unspaced, alternative settings) with its IT number;
           the kit must read point-group order x lattice points per cell equal to its own table's space-group order.
  cases:   primitive C-centred monoclinic cells with beta = 90 (two right angles, like the conventional cell), the
           operation lists spglib finds in that cell, and the atom count pymatgen's primitive structure holds.
"""
import json, pathlib
import numpy as np
import spglib
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

out = pathlib.Path(__file__).resolve().parents[3] / "Tests/DiffractoCoreTests/Fixtures/cif-symmetry-pins.json"

symbols = set()
for hall in range(1, 531):
    t = spglib.get_spacegroup_type(hall)
    for field in ("international", "international_full", "international_short"):
        for sym in getattr(t, field).split(" = "):
            symbols.add((t.number, sym))
symbols = sorted(symbols)


def xyz(rot, trans):
    parts = []
    for r, t in zip(rot, trans):
        terms = []
        for c, v in zip("xyz", r):
            if v == 1: terms.append(f"+{c}")
            elif v == -1: terms.append(f"-{c}")
            elif v != 0: terms.append(f"{v:+d}{c}")
        s = "".join(terms).lstrip("+")
        t = float(t) % 1.0
        if abs(t) > 1e-9:
            num = int(round(t * 12))
            from fractions import Fraction
            f = Fraction(num, 12)
            s += f"+{f.numerator}/{f.denominator}"
        parts.append(s)
    return ",".join(parts)


def cif(structure, symbol, number, ops, name, sites=None):
    p = structure.lattice.parameters
    lines = [f"data_{name}", f"_symmetry_space_group_name_H-M '{symbol}'"]
    if number is not None:
        lines.append(f"_symmetry_Int_Tables_number {number}")
    for tag, v in zip(["a", "b", "c"], p[:3]):
        lines.append(f"_cell_length_{tag} {v:.6f}")
    for tag, v in zip(["alpha", "beta", "gamma"], p[3:]):
        lines.append(f"_cell_angle_{tag} {v:.6f}")
    lines += ["loop_", "_symmetry_equiv_pos_as_xyz"] + [f"'{o}'" for o in ops]
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    for i, site in enumerate(structure if sites is None else [structure[k] for k in sites]):
        f = site.frac_coords
        lines.append(f"{site.specie.symbol}{i} {site.specie.symbol} {f[0]:.8f} {f[1]:.8f} {f[2]:.8f}")
    return "\n".join(lines) + "\n"


conv = Structure.from_spacegroup("C2/m", Lattice.monoclinic(5.0, 4.0, 3.0, 90.0), ["Si", "O"],
                                 [[0.1, 0, 0.2], [0.3, 0.2, 0.6]])
prim = SpacegroupAnalyzer(conv, symprec=1e-4).get_primitive_standard_structure()
cell = (prim.lattice.matrix, prim.frac_coords, [s.specie.Z for s in prim])
sym = spglib.get_symmetry(cell, symprec=1e-4)
ds = spglib.get_symmetry_dataset(cell, symprec=1e-4)
asym = sorted(set(int(k) for k in ds.equivalent_atoms))
prim_ops = sorted({xyz(r, t) for r, t in zip(sym["rotations"], sym["translations"])})
conv_ops_full = sorted({xyz(r, t) for r, t in zip(*[spglib.get_symmetry(
    (conv.lattice.matrix, conv.frac_coords, [s.specie.Z for s in conv]), 1e-4)[k] for k in ("rotations", "translations")])})
cases = [
    {"name": "primitive C2/m, beta 90, IT number", "text": cif(prim, "C 1 2/m 1", 12, prim_ops, "p1", asym),
     "atoms": len(prim), "asymmetricSites": len(asym), "operations": len(prim_ops)},
    {"name": "primitive C2/m, beta 90, symbol only", "text": cif(prim, "C 1 2/m 1", None, prim_ops, "p2", asym),
     "atoms": len(prim), "asymmetricSites": len(asym), "operations": len(prim_ops)},
]
json.dump({"pymatgen": "2026.9.24", "spglib": spglib.__version__, "symbols": symbols, "cases": cases,
           "primitiveParameters": list(prim.lattice.parameters), "primitiveOps": len(prim_ops),
           "conventionalOpsListed": len(conv_ops_full)}, open(out, "w"))
print(out, len(symbols), len(prim_ops), len(conv_ops_full), prim.lattice.parameters)
