"""orix pole-folding pins for the eleven Laue classes and both -3m settings (DiffractoKit 0.6.0), the reference for
`CrystalOrientationSymmetry.reduceDirection` and `sampleFundamentalZone`.

  ../tools/memcap.sh -- ../References/python/orix/bin/python -I -B tools/orix/direction-pins.py

Per class: orix's full Laue group as 3x3 matrices in the kit frame (F g F^T; see ipf-pins.py for F), poles in the kit frame
with `Vector3d.in_fundamental_sector` of each (also in the kit frame), and 1500 Fibonacci probe poles for the coverage check.
The kit's fundamental zone is not orix's (the kit takes the largest-(z,x,y) orbit member, or the app's native sector for
m-3m, 6/mmm and 2/m), so the Swift tests compare orbits: the kit's reduced pole must be in the orbit of the input and of orix's
folded pole under orix's own group, and every member of an orbit must reduce to the same pole.
Writes Tests/DiffractoOrientationTests/Fixtures/orix-direction-pins.json.
"""
import importlib.metadata as md
import json
from pathlib import Path
import numpy as np
from orix.quaternion import symmetry as S
from orix.vector import Vector3d

KIT = Path(__file__).resolve().parents[2]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261021
assert md.version('orix') == '0.15.0' and md.version('numpy') == '2.5.4'
CLASSES = [('-1', S.Ci), ('2/m', S.C2h), ('mmm', S.D2h), ('4/m', S.C4h), ('4/mmm', S.D4h), ('-3', S.S6),
           ('-3m', S.D3d), ('6/m', S.C6h), ('6/mmm', S.D6h), ('m-3', S.Th), ('m-3m', S.Oh), ('-3m1', S.D3d)]

def rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
FRAME = {'2/m': np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.0]]), '-3m': rz(-np.pi / 6)}

def fibonacci(n):
    i = np.arange(n)
    z = 1 - 2 * (i + .5) / n
    r = np.sqrt(1 - z * z)
    a = np.pi * (3 - np.sqrt(5)) * i
    return np.stack([r * np.cos(a), r * np.sin(a), z], axis=1)

def main():
    rng = np.random.default_rng(SEED)
    out = {'orix': md.version('orix'), 'seed': SEED, 'classes': []}
    for symbol, sym in CLASSES:
        F = FRAME.get(symbol, np.eye(3))
        group = np.array([(-1.0 if bool(sym[i].improper[0]) else 1.0) * sym[i].to_matrix().reshape(3, 3) for i in range(sym.size)])
        group = np.array([F @ g @ F.T for g in group])
        fs = sym.fundamental_sector
        poles = rng.normal(size=(300, 3))
        poles /= np.linalg.norm(poles, axis=1)[:, None]
        extra = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [0, 1, 0], [0, -1, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1], [0, 1, 1], [-1, 2, 3]]
        extra += [list(F @ v) for v in fs.vertices.data.reshape(-1, 3)] + [list(F @ fs.center.data.reshape(3))]
        poles = np.vstack([poles, np.array([np.array(e, float) / np.linalg.norm(e) for e in extra])])
        folded = np.asarray(Vector3d(poles @ F).in_fundamental_sector(sym.laue).data).reshape(-1, 3) @ F.T
        out['classes'].append({'symbol': symbol, 'group': group.tolist(), 'poles': poles.tolist(), 'orixFolded': folded.tolist(),
                               'probes': fibonacci(1500).tolist()})
    (FIX / 'orix-direction-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok', len(out['classes']))

if __name__ == '__main__':
    main()
