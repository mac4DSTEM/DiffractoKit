"""orix (Mis)orientation reduction pins for the eleven Laue classes (DiffractoKit 0.6.0).

  ../tools/memcap.sh -- ../References/DiffractoKit/orix-venv/bin/python -I -B tools/orix/disorientation-pins.py

Writes Tests/DiffractoOrientationTests/Fixtures/orix-disorientation-pins.json. Convention: the kit's matrix M (columns =
lab axes in crystal coordinates, py4DSTEM) is orix's Orientation matrix directly (orix applies the symmetry on the left:
reduced = S * O). orix works in its own settings; `frame` maps orix crystal coordinates to the kit's (see ipf-pins.py), so
M_orix = frame^T M_kit frame (the lab axes are relabelled with the crystal's, which keeps the angle to the
reference invariant) and M_kit = frame M_orix frame^T. Pairs use M_b M_a^T, which needs no lab relabelling. Singles: Orientation.reduce(); pairs: Orientation.angle_with and the
Misorientation((sym, sym)).reduce() route, which agree here and with brute force over the proper group (asserted).
"""
import importlib.metadata as md
import json
from pathlib import Path
import numpy as np
from orix.quaternion import Orientation, Misorientation, symmetry as S
from scipy.spatial.transform import Rotation

KIT = Path(__file__).resolve().parents[2]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261020
assert md.version('orix') == '0.15.0'
CLASSES = [('-1', S.Ci), ('2/m', S.C2h), ('mmm', S.D2h), ('4/m', S.C4h), ('4/mmm', S.D4h), ('-3', S.S6),
           ('-3m', S.D3d), ('6/m', S.C6h), ('6/mmm', S.D6h), ('m-3', S.Th), ('m-3m', S.Oh)]
def rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
FRAME = {'2/m': np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.0]]), '-3m': rz(-np.pi / 6)}
def orientation(M, sym):
    # Not Orientation.from_matrix: orix 0.15's Quaternion.from_matrix zeroes a small w near a half-turn (error 1.3e-5 at
    # |w| = 7e-6, found by the cross-check), so the quaternion comes from SciPy and orix's own to_matrix must give M back.
    q = Rotation.from_matrix(M).as_quat(scalar_first=True)
    for cand in (q, q * np.array([1, -1, -1, -1])):
        o = Orientation(cand, symmetry=sym)
        if np.abs(o.to_matrix().reshape(3, 3) - M).max() < 1e-12:
            return o
    raise AssertionError('quaternion convention')

def ang(R): return float(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))

def main():
    rng = np.random.default_rng(SEED)
    out = {'orix': md.version('orix'), 'seed': SEED, 'classes': []}
    for symbol, sym in CLASSES:
        F = FRAME.get(symbol, np.eye(3))
        G = sym.laue_proper_subgroup.to_matrix().reshape(-1, 3, 3)
        singles = list(Rotation.random(150, random_state=rng).as_matrix())
        singles += [np.eye(3)] + [Rotation.from_rotvec(np.array(a) * t).as_matrix() for a in ([0, 0, 1], [1, 1, 0], [1, 2, 3])
                                  for t in (0.01, 0.7, 3.0)]
        pairs = list(zip(Rotation.random(150, random_state=rng).as_matrix(), Rotation.random(150, random_state=rng).as_matrix()))
        rec = {'symbol': symbol, 'orixName': sym.name, 'singles': [], 'pairs': []}
        for Mk in singles:
            Mo = F.T @ Mk @ F   # lab axes relabelled too, so the angle to the reference is invariant
            red = orientation(Mo, sym).reduce()
            angles = sorted(ang(g @ Mo) for g in G)
            assert abs(float(red.angle[0]) - angles[0]) < 1e-9, symbol
            gap = angles[1] - angles[0] if len(angles) > 1 else 1.0
            rec['singles'].append({'matrix': Mk.tolist(), 'angle': float(red.angle[0]),
                                   'reduced': (F @ red.to_matrix().reshape(3, 3) @ F.T).tolist(), 'unique': bool(gap > 1e-6)})
        for Ma, Mb in pairs:
            oa = orientation(F.T @ Ma @ F, sym)
            ob = orientation(F.T @ Mb @ F, sym)
            angle = float(oa.angle_with(ob)[0])
            D = F.T @ Mb @ Ma.T @ F  # = (F^T Mb F)(F^T Ma F)^T
            brute = min(ang(g @ D) for g in G)
            m = Misorientation(orientation(D, sym).data); m.symmetry = (sym, sym)
            viaReduce = float(m.reduce().angle[0])
            assert abs(angle - brute) < 1e-8 and abs(viaReduce - brute) < 1e-8, (symbol, angle, brute, viaReduce)
            rec['pairs'].append({'a': Ma.tolist(), 'b': Mb.tolist(), 'angle': angle})
        out['classes'].append(rec)
    (FIX / 'orix-disorientation-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    print('ok')

if __name__ == '__main__':
    main()
