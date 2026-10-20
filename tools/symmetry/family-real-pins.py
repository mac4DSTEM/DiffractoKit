"""py4DSTEM `symmetry_reduce_directions` pins on REAL crystals (DiffractoKit 0.6.0): the operator list a real `Crystal` uses is
pymatgen's (`SpacegroupAnalyzer(structure).get_point_group_operations()`, duplicates removed), not the sorted list of the
harness in family-pins.py, and where several operators put a direction inside the plan sector the first in the list wins.

  ../tools/memcap.sh -- ../References/CuboKit/py4dstem-venv/bin/python -I -B tools/symmetry/family-real-pins.py

For each Laue class, representative space groups (centrosymmetric, a generic atom expanded with pymatgen's
`Structure.from_spacegroup`; hexagonal cells for the trigonal groups, R-3m and R-3c in hexagonal axes) build a real py4DSTEM
`Crystal` (pymatgen is imported first so `pymatgen_available` is true), `orientation_plan(zone_axis_range='auto')` chooses the
plan sector from the detected point group, and the pinned function runs on seeded random matrices. Records: the plan's
`symmetry_operators` (after duplicate removal and the -3 fibre doubling), the fractional rotation matrices pymatgen listed (the
order), the plan sector, and per matrix the families, the number of qualifying operators per column and boundary columns.
Needs pymatgen in the py4DSTEM venv (2025.1.9 there) and the spglib it brings. Writes
Tests/DiffractoOrientationTests/Fixtures/py4dstem-family-real-pins.json and Sources/DiffractoOrientation/PymatgenOperatorTable.swift
(the order pymatgen lists the point-group operations in, per class, from the first space group of each class above: P-1, P2/m, Pmmm,
P4/m, P4/mmm, P-3, P-31m, P-3m1, P6/m, P6/mmm, Pm-3, Pm-3m). Other space groups of a class that list the operations in another order
(origin choice 2, some centred groups) gave identical families on every pinned matrix.
"""
import contextlib
import importlib.metadata as md
import io
import json
from pathlib import Path
import subprocess
import warnings
import numpy as np
import pymatgen.core  # noqa: F401  (must precede the Crystal so that pymatgen_available is set)
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import py4DSTEM
from py4DSTEM.process.diffraction import Crystal
from py4DSTEM.process.diffraction.utils import Orientation
from scipy.spatial.transform import Rotation

warnings.filterwarnings('ignore')
KIT = Path(__file__).resolve().parents[2]
SUITE = KIT.parent
LOCK = 'f050d207851f8accd3aa550365baaf1f38230428'
clone = SUITE / 'References/CuboKit/py4DSTEM-f050d207'
assert subprocess.check_output(['git', '-C', str(clone), 'rev-parse', 'HEAD'], text=True).strip() == LOCK
assert py4DSTEM.__version__ == '0.14.19'
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261021
# class -> space groups. The first of each class is the one whose order the kit ships ("primary"); the others show which orders
# differ (the origin-choice-2 and some centred groups list the operations in another order).
GROUPS = {'-1': [2], '2/m': [10, 11, 12, 13, 14, 15], 'mmm': [47, 51, 62, 63, 65, 69, 71, 74, 70],
          '4/m': [83, 84, 87, 85, 86, 88], '4/mmm': [123, 127, 131, 136, 139, 140, 129, 141, 142], '-3': [147, 148],
          '-3m': [162, 163], '-3m1': [164, 165, 166, 167], '6/m': [175, 176], '6/mmm': [191, 192, 193, 194],
          'm-3': [200, 202, 204, 205, 206, 201, 203], 'm-3m': [221, 225, 229, 230, 227]}

def lattice(sg):
    if sg < 3: return Lattice.from_parameters(3.1, 3.7, 4.3, 80, 95, 100)
    if sg < 16: return Lattice.monoclinic(3.1, 3.7, 4.3, 100)
    if sg < 75: return Lattice.orthorhombic(3.1, 3.7, 4.3)
    if sg < 143: return Lattice.tetragonal(3.1, 4.3)
    if sg < 195: return Lattice.hexagonal(3.1, 4.3)
    return Lattice.cubic(3.5)

def main():
    rng = np.random.default_rng(SEED)
    mats = list(Rotation.random(120, random_state=rng).as_matrix())
    mats += [np.eye(3)] + [Rotation.from_rotvec(np.array(a) * t).as_matrix()
                           for a in ([0, 0, 1], [1, 1, 0], [1, 1, 1], [1, 2, 3], [1, 0, 0], [0, 1, 0]) for t in (0.3, 1.1, 2.0, 3.1, np.pi)]
    out = {'lock': LOCK, 'py4DSTEM': py4DSTEM.__version__, 'pymatgen': md.version('pymatgen'), 'spglib': md.version('spglib'),
           'numpy': np.__version__, 'seed': SEED, 'cases': []}
    for key, numbers in GROUPS.items():
        for sg in numbers:
            s = Structure.from_spacegroup(sg, lattice(sg), ['Si'], [[0.1234, 0.2345, 0.3456]])
            cell = [s.lattice.a, s.lattice.b, s.lattice.c, s.lattice.alpha, s.lattice.beta, s.lattice.gamma]
            cr = Crystal(np.array([x.frac_coords for x in s]), np.array([14] * len(s)), cell)
            assert cr.pymatgen_available
            with contextlib.redirect_stdout(io.StringIO()):
                cr.orientation_plan(zone_axis_range='auto', angle_step_zone_axis=45, angle_step_in_plane=90,
                                    calculate_correlation_array=False, progress_bar=False)
            sga = cr.pointgroup
            symbol = sga.get_point_group_symbol()
            assert symbol == key.replace('-3m1', '-3m'), (sg, symbol)
            order, seen = [], set()
            for o in sga.get_point_group_operations():
                m = tuple(int(x) for x in np.rint(o.rotation_matrix).flatten())
                if m not in seen:
                    seen.add(m); order.append(list(m))
            cases = []
            for M in mats:
                o = Orientation(1); o.matrix[0] = M
                o = cr.symmetry_reduce_directions(o)
                hits, boundary = [], []
                for a in range(3):
                    coords = np.sum(cr.symmetry_reduction * M[:, a][None, :, None], axis=1)
                    hits.append(int(np.count_nonzero(np.all(coords >= 0, axis=1))))
                    boundary.append(bool(np.any((coords.min(axis=1) >= -1e-12) & (coords.min(axis=1) <= 1e-12))))
                cases.append({'matrix': M.tolist(), 'family': o.family[0].tolist(), 'hits': hits, 'boundary': boundary})
            lat = cr.lat_real
            out['cases'].append({'class': key, 'spaceGroup': sg, 'symbol': sga.get_space_group_symbol(), 'primary': sg == numbers[0],
                'lattice': {'a': cell[0], 'b': cell[1], 'c': cell[2], 'alpha': cell[3], 'beta': cell[4], 'gamma': cell[5]},
                'sites': [{'z': 14, 'abc': [float(t) for t in x.frac_coords]} for x in s],
                'orderFractional': order, 'operators': cr.symmetry_operators.tolist(),
                'range': cr.orientation_zone_axis_range.tolist(), 'cases': cases})
            print(key, sg, len(order), len(cr.symmetry_operators), flush=True)
    first = {}
    for c in out['cases']:
        first.setdefault(c['class'], c['orderFractional'])
    case = {'-1': ('oneBar', None), '2/m': ('twoOverM', None), 'mmm': ('mmm', None), '4/m': ('fourOverM', None),
            '4/mmm': ('fourOverMmm', None), '-3': ('threeBar', None),
            '-3m1': ('threeBarM', 'threeBarMOne'), '-3m': ('threeBarM', None),
            # (-3m1 first: the pattern after it is the default -31m)
            '6/m': ('sixOverM', None), '6/mmm': ('sixOverMmm', None),
            'm-3': ('mThreeBar', None), 'm-3m': ('mThreeBarM', None)}
    lines = ['// Generated by tools/symmetry/family-real-pins.py. Do not edit.',
             '// The order in which pymatgen (SpacegroupAnalyzer.get_point_group_operations, as py4DSTEM 0.14.19 reads it) lists the',
             '// operations of the point group of a centrosymmetric structure: rotation matrices in the fractional basis, row major.',
             'import DiffractoCore', '',
             'enum PymatgenOperatorTable {',
             '    static func fractional(_ laueClass: LaueClass, _ setting: ThreeBarMSetting?) -> [[Int]] {',
             '        switch (laueClass, setting) {']
    for key, (cls, setting) in case.items():
        pattern = '(.%s, .%s)' % (cls, setting) if setting else '(.%s, _)' % cls
        lines.append('        case %s:' % pattern)
        lines.append('            return [' + ', '.join('[' + ','.join(str(x) for x in m) + ']' for m in first[key]) + ']')
    lines += ['        }', '    }', '}', '']
    (KIT / 'Sources/DiffractoOrientation/PymatgenOperatorTable.swift').write_text('\n'.join(lines))
    (FIX / 'py4dstem-family-real-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')

if __name__ == '__main__':
    main()
