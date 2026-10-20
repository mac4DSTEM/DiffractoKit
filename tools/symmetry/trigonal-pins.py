"""Pins for the two -3m settings (DiffractoKit 0.6.0): the point-group operations pymatgen/spglib give for real trigonal
structures, in the order they give them, in hexagonal fractional coordinates and in the kit's Cartesian frame.

  ../tools/memcap.sh -- ../References/DiffractoKit/pymatgen-venv/bin/python -I -B tools/symmetry/trigonal-pins.py

For space groups 162, 163 (-31m type) and 164, 165, 166, 167 (-3m1 type, 166 and 167 in hexagonal axes) a generic-position atom
is expanded with `Structure.from_spacegroup`, the lattice is hexagonal (a 3.1, c 4.3), `SpacegroupAnalyzer.get_point_group_operations()`
lists the operations (this is the list py4DSTEM's `Crystal` reads), and the twofold axes are classified from the Cartesian matrices
(angle of the axis from x: 0, 60, 120 degrees is -3m1, 30, 90, 150 degrees is -31m). Writes
Tests/DiffractoOrientationTests/Fixtures/trigonal-pins.json.
"""
import importlib.metadata as md
import json
from pathlib import Path
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

KIT = Path(__file__).resolve().parents[2]
OUT = KIT / 'Tests/DiffractoOrientationTests/Fixtures/trigonal-pins.json'
assert md.version('pymatgen') == '2026.9.24' and md.version('numpy') == '2.5.4'
H = np.array([[1, 0, 0], [-.5, np.sqrt(3) / 2, 0], [0, 0, 1]])   # rows: a, b, c in the kit frame

def twofold_angle_deg(R):
    w, v = np.linalg.eig(R)
    axis = np.real(v[:, np.argmin(np.abs(w - 1))])
    if abs(axis[2]) > 1e-9:
        return None
    return float(np.rad2deg(np.arctan2(axis[1], axis[0])) % 180.0)

def main():
    out = {'pymatgen': md.version('pymatgen'), 'spglib': md.version('spglib'), 'cases': []}
    for number in (162, 163, 164, 165, 166, 167):
        structure = Structure.from_spacegroup(number, Lattice.hexagonal(3.1, 4.3), ['Si'], [[0.1234, 0.2345, 0.3456]])
        sga = SpacegroupAnalyzer(structure, symprec=1e-3)
        assert sga.get_point_group_symbol() == '-3m'
        ops = sga.get_point_group_operations()
        order, seen = [], set()
        for o in ops:
            m = np.rint(o.rotation_matrix).astype(int)
            assert np.allclose(m, o.rotation_matrix)
            key = tuple(int(x) for x in m.flatten())
            if key not in seen:
                seen.add(key); order.append(m)
        assert len(order) == 12, (number, len(order))
        cart = [H.T @ m @ np.linalg.inv(H.T) for m in order]   # same transform as tools/symmetry/pins.py
        proper = [(m, c) for m, c in zip(order, cart) if np.linalg.det(m) > 0]
        angles = sorted({round(a, 6) for m, c in proper if np.trace(c) < 0 and (a := twofold_angle_deg(c)) is not None})
        setting = '-3m1' if angles == [0.0, 60.0, 120.0] else '-31m' if angles == [30.0, 90.0, 150.0] else None
        assert setting is not None, (number, angles)
        out['cases'].append({
            'spaceGroup': number, 'symbol': sga.get_space_group_symbol(), 'setting': setting, 'twofoldAnglesDeg': angles,
            'lattice': {'a': 3.1, 'b': 3.1, 'c': 4.3, 'alpha': 90.0, 'beta': 90.0, 'gamma': 120.0},
            'sites': [{'z': int(s.specie.Z), 'abc': [float(x) for x in s.frac_coords]} for s in structure],
            'operatorsFractional': [m.tolist() for m in order],               # all 12, pymatgen's order
            'operatorsCartesian': [c.tolist() for c in cart]})
    OUT.write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok', [(c['spaceGroup'], c['setting'], len(c['sites'])) for c in out['cases']])

if __name__ == '__main__':
    main()
