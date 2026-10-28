"""orix Misorientation / Orientation pins (DiffractoKit 0.10.0).

  ../References/python/orix/bin/python -I -B tools/orix/misorientation-pins.py

Writes Tests/DiffractoOrientationTests/Fixtures/orix-misorientation-pins.json. The oracle is orix 0.15.0:
  * `Misorientation(M, symmetry=(Gl, Gr))`: `equivalent(grain_exchange)` = Gr * M * Gl, `reduce()` (the smallest rotation angle).
    The kit's `Misorientation(matrix, source: Gl, target: Gr)` is the same object: the matrix maps source-crystal vectors to
    target-crystal vectors, the equivalents are g_target R g_source. Asserted here: the matrix of each orix equivalent equals the
    product of the matrices of its factors.
  * `Orientation.in_euler_fundamental_region()` for the eleven Laue classes (and -3m in both settings).
The symmetry groups are built from the kit's own generators (CrystalOrientationSymmetry.init, copied below), so orix and the kit
share the settings: unique axis b for 2/m, hexagonal a = x, b at 120 degrees, -31m (twofold perpendicular to a) by default.
Equivalents are compared as sets (atol 1e-9); the reduced matrix only where the smallest angle is unique (gap above 1e-6);
angles at atol 1e-4 (orix's own test tolerance).
The beta'' case (task 35): R0 from tools/beta-variants/planted.py (numpy-only derivation, copied below), source = 2/m (b unique),
target = m-3m: orix's equivalent() has 24 distinct elements and 12 classes under the source group; the reduced angle is arctan(2/3).
"""
import importlib.metadata as md
import itertools
import json
from pathlib import Path
import numpy as np
from orix.quaternion import Misorientation, Orientation, Rotation, Symmetry, symmetry as S
from scipy.spatial.transform import Rotation as SR

KIT = Path(__file__).resolve().parents[2]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261022
assert md.version('orix') == '0.15.0'

def im(*rows): return np.array(rows, dtype=float)
X = im([1, 0, 0], [0, -1, 0], [0, 0, -1]); Y = im([-1, 0, 0], [0, 1, 0], [0, 0, -1])
FOUR = im([0, -1, 0], [1, 0, 0], [0, 0, 1]); THREE = im([0, -1, 0], [1, -1, 0], [0, 0, 1])
SIX = im([1, -1, 0], [1, 0, 0], [0, 0, 1]); DIAG = im([0, 0, 1], [1, 0, 0], [0, 1, 0])
TWO31 = im([0, -1, 0], [-1, 0, 0], [0, 0, -1]); TWO3M1 = im([1, -1, 0], [0, -1, 0], [0, 0, -1])
GEN = {'-1': [], '2/m': [Y], 'mmm': [X, Y], '4/m': [FOUR], '4/mmm': [FOUR, X], '-3': [THREE], '-31m': [THREE, TWO31],
       '-3m1': [THREE, TWO3M1], '6/m': [SIX], '6/mmm': [SIX, TWO3M1], 'm-3': [X, DIAG], 'm-3m': [FOUR, DIAG]}
HEX = {'-3', '-31m', '-3m1', '6/m', '6/mmm'}
HEXFRAME = np.array([[1, -0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1.0]])

def closure(gens):
    ops = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for a in list(ops):
            for g in gens:
                p = g @ a
                if not any(np.allclose(p, o) for o in ops):
                    ops.append(p); changed = True
    return ops

def cartesian(symbol):
    f = HEXFRAME if symbol in HEX else np.eye(3)
    return [f @ r @ np.linalg.inv(f) for r in closure(GEN[symbol])]

def rot(M):
    q = SR.from_matrix(M).as_quat(scalar_first=True)
    for cand in (q, q * np.array([1, -1, -1, -1])):
        r = Rotation(cand)
        if np.abs(r.to_matrix().reshape(3, 3) - M).max() < 1e-12:
            return r
    raise AssertionError('quaternion convention')

def symmetry(symbol):
    """The orix symmetry (proper operators of the kit's setting plus inversion, so orix names the proper subgroup) and the proper matrices."""
    ops = cartesian(symbol)
    if symbol == '-1':
        return S.Ci, ops
    q = np.concatenate([rot(m).data for m in ops])
    sym = Symmetry(np.concatenate([q, q]))
    sym.improper = np.array([0] * len(ops) + [1] * len(ops))
    return sym, ops

def ang(R): return float(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))

def mats(r): return r.to_matrix().reshape(-1, 3, 3)

def main():
    rng = np.random.default_rng(SEED)
    out = {'orix': md.version('orix'), 'seed': SEED, 'pairs': [], 'euler': [], 'beta': None}
    # NOTE: orix symmetries here are the PROPER operators only: the kit's Misorientation works on rotations; inversion changes nothing.
    pairs = [('m-3m', '2/m'), ('2/m', 'm-3m'), ('m-3m', '6/mmm'), ('6/mmm', '2/m'), ('4/mmm', 'mmm'), ('m-3', '-3'), ('-31m', '-3m1'),
             ('-1', 'm-3m'), ('mmm', 'mmm'), ('6/m', '4/m'), ('m-3m', 'm-3m'), ('-3m1', '6/mmm')]
    for a, b in pairs:
        sa, oa = symmetry(a)
        sb, ob = symmetry(b)
        gl = sa.proper_subgroup
        gr = sb.proper_subgroup
        la = mats(gl); lb = mats(gr)
        rec = {'source': a, 'target': b, 'sourceOrder': len(la), 'targetOrder': len(lb), 'items': []}
        Ms = list(SR.random(40, random_state=rng).as_matrix()) + [np.eye(3)]
        keepSets = 2          # only the first few items carry the full list of equivalents (file size)
        for index, M in enumerate(Ms):
            m = Misorientation(rot(M).data); m.symmetry = (gl, gr)
            eq = m.equivalent()
            em = mats(eq)
            brute = [lb_i @ M @ la_j for lb_i in lb for la_j in la]
            # orix's product order equals the matrix product order
            assert all(any(np.allclose(e, bm, atol=1e-9) for bm in brute) for e in em)
            uniq = []
            for e in em:
                if not any(np.allclose(e, u, atol=1e-9) for u in uniq):
                    uniq.append(e)
            red = m.reduce()
            angles = sorted(ang(g) for g in brute)
            gap = [x for x in np.diff(angles) if x > 1e-9]
            unique = bool(len(angles) < 2 or (angles[1] - angles[0]) > 1e-6)
            orix_min = abs(float(red.angle[0]) - angles[0]) < 1e-8
            item = {'matrix': M.tolist(), 'equivalentCount': len(uniq), 'angle': float(red.angle[0]), 'minAngle': angles[0], 'orixIsMinimal': bool(orix_min),
                    'unique': unique, 'reduced': mats(red)[0].tolist()}
            if index < keepSets or index == len(Ms) - 1:
                item['equivalents'] = [u.tolist() for u in uniq]
            if a == b:
                ge = m.equivalent(grain_exchange=True)
                guniq = []
                for e in mats(ge):
                    if not any(np.allclose(e, u, atol=1e-9) for u in guniq):
                        guniq.append(e)
                item['grainExchangeCount'] = len(guniq)
            rec['items'].append(item)
        out['pairs'].append(rec)

    # in_euler_fundamental_region
    for symbol in ['-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-31m', '-3m1', '6/m', '6/mmm', 'm-3', 'm-3m']:
        s, ops = symmetry(symbol)
        pg = s.proper_subgroup
        rec = {'symbol': symbol, 'orixProperName': pg.name, 'region': list(pg.euler_fundamental_region), 'items': []}
        for M in list(SR.random(60, random_state=rng).as_matrix()) + [np.eye(3)]:
            o = Orientation(rot(M).data, symmetry=s)
            eu = o.in_euler_fundamental_region()
            rec['items'].append({'matrix': M.tolist(), 'euler': [float(v) for v in eu[0]]})
        out['euler'].append(rec)

    # beta'' in Al
    s13 = np.sqrt(13.0)
    ea = np.array([3.0, 2.0, 0.0]) / s13; eb = np.array([0.0, 0.0, 1.0])
    R0 = np.column_stack([ea, eb, np.cross(ea, eb)])
    sp, _ = symmetry('2/m'); sm, _ = symmetry('m-3m')
    gp = sp.proper_subgroup; gm = sm.proper_subgroup
    m = Misorientation(rot(R0).data); m.symmetry = (gp, gm)
    eq = mats(m.equivalent())
    uniq = []
    for e in eq:
        if not any(np.allclose(e, u, atol=1e-9) for u in uniq):
            uniq.append(e)
    red = m.reduce()
    out['beta'] = {'R0': R0.tolist(), 'equivalentCount': len(uniq), 'reducedAngle': float(red.angle[0]),
                   'reducedAngleClosedForm': float(np.arctan2(2, 3))}
    (FIX / 'orix-misorientation-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('orix reduce() is not the minimum angle in', {f"{p['source']}>{p['target']}": sum(not i['orixIsMinimal'] for i in p['items']) for p in out['pairs']})
    print('ok', [(p['source'], p['target'], p['items'][0]['equivalentCount']) for p in out['pairs']],
          [(e['symbol'], e['orixProperName'], e['region']) for e in out['euler']], out['beta']['equivalentCount'], out['beta']['reducedAngle'])

if __name__ == '__main__':
    main()
