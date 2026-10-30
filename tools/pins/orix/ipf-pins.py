"""orix IPF colour key and fundamental-sector pins for the eleven Laue classes (DiffractoKit 0.6.0).

Run from the kit root with the pinned orix venv (References/python/orix):
  ../References/python/orix/bin/python -I -B tools/pins/orix/ipf-pins.py

Writes Tests/DiffractoOrientationTests/Fixtures/orix-ipf-pins.json and Sources/DiffractoOrientation/OrixIPFSectors.swift.
orix is the oracle (MTEX is MATLAB and is not used). orix works in its own settings (2/m with the twofold on z, -3m with
twofolds on x); `frame` is the rotation that maps orix coordinates into the kit's documented standard setting:
  2/m   kit (x,y,z) = (Y,Z,X)_orix (cyclic: the kit's unique axis b = y is orix's z)
  -3m   Rz(-30 deg) (the kit's pymatgen setting has twofolds at -30, 30, 90 degrees; orix's at 0, 60, 120)
  -3m1  identity: the kit's -3m1 setting (P-3m1, R-3m; twofolds on a) is orix's D3d itself. It is the last entry, so the other
        classes' seeded draws are unchanged. Its kit operators come from tools/pins/orientation-symmetry/trigonal-pins.py (pymatgen).
  others identity.
The generator asserts that frame * orix proper operators == the kit's operators as sets (symmetry-pins.json).
"""
import importlib.metadata as md
import json
from pathlib import Path
import sys
import numpy as np
from orix.quaternion import symmetry as S
from orix.plot.direction_color_keys import DirectionColorKeyTSL
from orix.vector import Vector3d

KIT = Path(__file__).resolve().parents[3]
FIX = KIT / 'Tests/DiffractoOrientationTests/Fixtures'
SEED = 20261020
assert md.version('orix') == '0.15.0'
assert md.version('numpy') == '2.5.4'

CLASSES = [('-1', S.Ci), ('2/m', S.C2h), ('mmm', S.D2h), ('4/m', S.C4h), ('4/mmm', S.D4h), ('-3', S.S6),
           ('-3m', S.D3d), ('6/m', S.C6h), ('6/mmm', S.D6h), ('m-3', S.Th), ('m-3m', S.Oh),
           ('-3m1', S.D3d)]

def rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

FRAME = {'2/m': np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.0]]), '-3m': rz(-np.pi / 6)}

def same_set(a, b):
    return len(a) == len(b) and all(any(np.allclose(x, y, atol=1e-12) for y in b) for x in a)

def main():
    kit = {g['symbol']: np.array(g['operators']) for g in
           json.load(open(FIX / 'symmetry-pins.json'))['groups']}
    tri = json.load(open(FIX / 'trigonal-pins.json'))['cases']
    kit['-3m1'] = np.array([c for c in next(t for t in tri if t['spaceGroup'] == 164)['operatorsCartesian']
                            if np.linalg.det(c) > 0])
    rng = np.random.default_rng(SEED)
    out = {'orix': md.version('orix'), 'numpy': md.version('numpy'), 'seed': SEED, 'classes': []}
    swift = []
    for symbol, sym in CLASSES:
        proper = sym.laue_proper_subgroup.to_matrix().reshape(-1, 3, 3)
        F = FRAME.get(symbol, np.eye(3))
        assert same_set([F @ m @ F.T for m in proper], list(kit[symbol])), symbol
        fs = sym.fundamental_sector
        key = DirectionColorKeyTSL(sym)
        normals = fs.data.reshape(-1, 3)
        center = fs.center.data.reshape(3)
        vertices = fs.vertices.data.reshape(-1, 3)
        # Kit-frame inputs: random signed unit vectors, the axes, the sector vertices and centre, edge midpoints.
        d = rng.normal(size=(240, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        extra = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]
        extra += [list(F @ v) for v in vertices] + [list(F @ center)]
        # Interior points near the vertices and edges (boundary-adjacent but not on a tie).
        for v in vertices:
            extra.append(list(F @ (v + 1e-4 * normals.mean(axis=0) + 1e-4 * center)))
        dk = np.vstack([d, np.array(extra, dtype=float)])
        do = dk @ F          # rows: F^T d  (orix coordinates)
        v = Vector3d(do)
        folded = np.asarray(v.in_fundamental_sector(sym.laue).data).reshape(-1, 3)
        rgb = key.direction2color(v).reshape(-1, 3)
        out['classes'].append({'symbol': symbol, 'orixName': sym.name, 'frame': F.tolist(), 'normals': normals.tolist(),
            'center': center.tolist(), 'vertices': vertices.tolist(), 'directions': dk.tolist(),
            'folded': (folded @ F.T).tolist(), 'rgb': rgb.tolist()})
        group = np.array([(-1.0 if bool(sym[i].improper[0]) else 1.0) * sym[i].to_matrix().reshape(3, 3)
                          for i in range(sym.size)])
        if symbol != '-3m1':   # the sector table is shared with -3m; only the frame differs
            swift.append((symbol, F, normals, center, vertices, group))
    (FIX / 'orix-ipf-pins.json').write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    case = {'-1': 'oneBar', '2/m': 'twoOverM', 'mmm': 'mmm', '4/m': 'fourOverM', '4/mmm': 'fourOverMmm', '-3': 'threeBar',
            '-3m': 'threeBarM', '6/m': 'sixOverM', '6/mmm': 'sixOverMmm', 'm-3': 'mThreeBar', 'm-3m': 'mThreeBarM'}
    def vec(x): return 'SIMD3(%s)' % ','.join(repr(float(t)) for t in x)
    def mat(F): return 'simd_double3x3(rows:[%s])' % ','.join(vec(r) for r in F)
    lines = ['// Generated by tools/pins/orix/ipf-pins.py from orix %s (fundamental sectors, taken from MTEX). Do not edit.' % md.version('orix'),
             'import simd', 'import DiffractoCore', '',
             '/// One Laue class\'s orix IPF fundamental sector, in orix\'s own frame, with the rotation `frame` that maps orix',
             '/// coordinates into the kit\'s standard setting.',
             'struct OrixIPFSector {',
             '    let frame: simd_double3x3', '    let normals: [SIMD3<Double>]', '    let center: SIMD3<Double>',
             '    let vertices: [SIMD3<Double>]',
             '    /// The full Laue group in orix\'s element order (proper elements first); ties pick the first.',
             '    let group: [simd_double3x3]', '',
             '    static func sector(for laueClass: LaueClass) -> OrixIPFSector {',
             '        switch laueClass {']
    for symbol, F, n, c, vs, grp in swift:
        lines += ['        case .%s:' % case[symbol],
                  '            return OrixIPFSector(frame: %s,' % mat(F),
                  '                normals: [%s],' % ','.join(vec(x) for x in n),
                  '                center: %s,' % vec(c),
                  '                vertices: [%s],' % ','.join(vec(x) for x in vs),
                  '                group: [%s])' % ','.join(mat(x) for x in grp)]
    lines += ['        }', '    }', '}', '']
    (KIT / 'Sources/DiffractoOrientation/OrixIPFSectors.swift').write_text('\n'.join(lines))
    print('ok', len(out['classes']))

if __name__ == '__main__':
    main()
