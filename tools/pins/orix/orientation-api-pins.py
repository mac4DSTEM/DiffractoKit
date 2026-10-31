"""orix Orientation API pins: `from_euler(symmetry=, direction=)`, `dot`, `dot_outer`, `get_distance_matrix`.

  ../References/python/orix/bin/python -I -B tools/pins/orix/orientation-api-pins.py

Writes Tests/DiffractoOrientationTests/Fixtures/orix-orientation-api-pins.json. The oracle is orix 0.15.0. The symmetry groups
and the quaternion helper come from misorientation-pins.py next to this file (the kit's own generators and settings).
  * `Orientation.from_euler(euler, symmetry, direction)`: the rotation matrix, Bunge angles in radians, both directions,
    random triplets plus the gimbal-lock and half-turn edges.
  * For each of the twelve class/setting entries: six orientations A and five B (A[1] and B[0] symmetry equivalents of A[0],
    the last operator, so the dot is 1 and the distance 0); `A[:5].dot(B)`, `A.dot_outer(B)` and `A.get_distance_matrix()`.
Tolerance: atol 1e-4 for from_euler (test_from_euler_symmetry), np.allclose defaults (rtol 1e-5, atol 1e-8) for the rest
(orix test_orientation.py).
"""
import importlib.metadata as md
import importlib.util
import json
from pathlib import Path
import numpy as np
from orix.quaternion import Orientation
from scipy.spatial.transform import Rotation as SR

HERE = Path(__file__).resolve().parent
KIT = Path(__file__).resolve().parents[3]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261031
assert md.version('orix') == '0.15.0'

spec = importlib.util.spec_from_file_location('misorientation_pins', HERE / 'misorientation-pins.py')
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)

SYMBOLS = ['-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-31m', '-3m1', '6/m', '6/mmm', 'm-3', 'm-3m']


def main():
    rng = np.random.default_rng(SEED)
    out = {'orix': md.version('orix'), 'seed': SEED, 'euler': [], 'sets': []}

    edges = [[0, 0, 0], [0.3, 0, 0.4], [0.3, np.pi, 0.4], [np.pi, np.pi / 2, np.pi], [2 * np.pi - 1e-3, 1e-3, 1e-3], [1.0, 2.0, 5.5]]
    triplets = edges + (rng.random((30, 3)) * np.array([2 * np.pi, np.pi, 2 * np.pi])).tolist()
    s_cubic, _ = P.symmetry('m-3m')
    for eu in triplets:
        rec = {'euler': [float(v) for v in eu]}
        for direction in ('lab2crystal', 'crystal2lab'):
            o = Orientation.from_euler(np.array(eu), symmetry=s_cubic, direction=direction)
            assert o.symmetry.name == s_cubic.name
            rec[direction] = P.mats(o)[0].tolist()
        out['euler'].append(rec)

    for symbol in SYMBOLS:
        s, ops = P.symmetry(symbol)
        A = list(SR.random(6, random_state=rng).as_matrix())
        A[1] = ops[-1] @ A[0]           # an equivalent of A[0]: distance 0 inside the distance matrix
        B = [ops[len(ops) // 2] @ A[0]] + list(SR.random(4, random_state=rng).as_matrix())
        oa = Orientation(np.concatenate([P.rot(m).data for m in A]), symmetry=s)
        ob = Orientation(np.concatenate([P.rot(m).data for m in B]), symmetry=s)
        dot = oa[:5].dot(ob)
        outer = oa.dot_outer(ob)
        dist = oa.get_distance_matrix()
        assert outer.shape == (6, 5) and dist.shape == (6, 6)
        assert np.allclose(dot, np.diag(outer[:5]))
        assert float(dist[0, 0]) == 0 and float(dist[0, 1]) == 0 and abs(float(outer[0, 0]) - 1) < 1e-12 and abs(float(oa[0].angle_with(ob[0])[0])) < 1e-6
        out['sets'].append({'symbol': symbol, 'A': [m.tolist() for m in A], 'B': [m.tolist() for m in B],
                            'dot': [float(v) for v in dot], 'dotOuter': outer.tolist(), 'distanceMatrix': dist.tolist()})

    (FIX / 'orix-orientation-api-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok', len(out['euler']), 'Euler triplets,', len(out['sets']), 'symmetry sets')


if __name__ == '__main__':
    main()
