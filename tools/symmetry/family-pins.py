"""py4DSTEM `symmetry_reduce_directions` pins (DiffractoKit 0.6.0): the three direction families of an orientation matrix.

  ../tools/memcap.sh -- ../References/CuboKit/py4dstem-venv/bin/python -I -B tools/symmetry/family-pins.py

Needs symmetry-pins.json (tools/symmetry/pins.py) for the proper operators (kit frame, sorted by their fractional
matrices) and the plan sector `range` of each class. Uses the same harness as pins-py4dstem.py: py4DSTEM's own
orientation_plan builds the zone-axis range (and the 180 degree fibre flip for -3), symmetry_operators = ops then -ops,
symmetry_reduction as py4DSTEM computes it, then the pinned function runs unchanged. Writes
Tests/DiffractoOrientationTests/Fixtures/py4dstem-family-pins.json. Also records, per matrix column, how many operators
satisfy the in-range test (the port picks the first) and which columns lie on a sector boundary within 1e-12.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import types
import numpy as np
import py4DSTEM
from py4DSTEM.process.diffraction import Crystal
from py4DSTEM.process.diffraction.crystal_ACOM import orientation_ranges
from py4DSTEM.process.diffraction.utils import Orientation
from scipy.spatial.transform import Rotation

KIT = Path(__file__).resolve().parents[2]
SUITE = KIT.parent
LOCK = 'f050d207851f8accd3aa550365baaf1f38230428'
clone = SUITE / 'References/CuboKit/py4DSTEM-f050d207'
assert subprocess.check_output(['git', '-C', str(clone), 'rev-parse', 'HEAD'], text=True).strip() == LOCK
assert py4DSTEM.__version__ == '0.14.19'
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261020

def main():
    pins = json.load(open(FIX / 'symmetry-pins.json'))
    rng = np.random.default_rng(SEED)
    out = {'lock': LOCK, 'version': py4DSTEM.__version__, 'numpy': np.__version__, 'seed': SEED, 'classes': []}
    for group in pins['groups']:
        symbol = group['symbol']
        zone, axis, angles = orientation_ranges[symbol]
        hexagonal = symbol in ['-3', '-3m', '6/m', '6/mmm']
        cr = Crystal([[0, 0, 0]], [14], [1, 1, 1, 90, 90, 120] if hexagonal else [1, 1, 1, 90, 90, 90])
        cr.pymatgen_available = False
        with contextlib.redirect_stdout(io.StringIO()):
            cr.orientation_plan(zone_axis_range=np.array(zone) if isinstance(zone, list) else zone,
                fiber_axis=axis, fiber_angles=angles, angle_step_zone_axis=45, angle_step_in_plane=90,
                calculate_correlation_array=False, progress_bar=False)
        ops = np.array(group['operators'])
        cr.symmetry_operators = np.concatenate((ops, -ops))
        inv = np.linalg.inv(cr.orientation_zone_axis_range)
        cr.symmetry_reduction = np.array([(inv.T @ op).T for op in cr.symmetry_operators])
        flip = bool(cr.orientation_fiber_angles is not None and abs(cr.orientation_fiber_angles[0] - 180) < 1e-3)
        flipInv = None
        if flip:
            flipped = cr.orientation_zone_axis_range.copy(); flipped[0] *= -1
            flipInv = np.linalg.inv(flipped)
            cr.symmetry_reduction = np.concatenate((cr.symmetry_reduction,
                np.array([(flipInv.T @ op).T for op in cr.symmetry_operators])))
            cr.symmetry_operators = np.tile(cr.symmetry_operators, (2, 1, 1))
        mats = list(Rotation.random(200, random_state=rng).as_matrix())
        mats += [np.eye(3)] + [Rotation.from_rotvec(np.array(a) * t).as_matrix()
                               for a in ([0, 0, 1], [1, 1, 0], [1, 1, 1], [1, 2, 3]) for t in (0.3, 1.1, 2.0, 3.1)]
        cases = []
        for M in mats:
            o = Orientation(1); o.matrix[0] = M
            o = cr.symmetry_reduce_directions(o)
            hits = [int(np.count_nonzero(np.all(np.sum(cr.symmetry_reduction * M[:, a][None, :, None], axis=1) >= 0, axis=1)))
                    for a in range(3)]
            # A column is a boundary column when some operator puts it within 1e-12 of a sector wall (all its coordinates
            # >= -1e-12 and the smallest within 1e-12 of zero): whether that operator passes the exact `>= 0` test is
            # floating-point noise of NumPy's BLAS product, so such columns are only checked loosely.
            boundary = []
            for a in range(3):
                coords = np.sum(cr.symmetry_reduction * M[:, a][None, :, None], axis=1)
                boundary.append(bool(np.any((coords.min(axis=1) >= -1e-12) & (coords.min(axis=1) <= 1e-12))))
            cases.append({'matrix': M.tolist(), 'family': o.family[0].tolist(), 'hits': hits, 'boundary': boundary})
        out['classes'].append({'symbol': symbol, 'range': cr.orientation_zone_axis_range.tolist(), 'flip': flip,
            'basis': inv.T.tolist(), 'flipBasis': flipInv.T.tolist() if flip else None,
            'operatorCount': len(cr.symmetry_operators), 'cases': cases})
    (FIX / 'py4dstem-family-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok')

if __name__ == '__main__':
    main()
