"""Generate symmetry fixtures using the pinned pymatgen environment, under a memory cap.

From this worktree:
  ../../References/python/pymatgen/bin/python -I -B tools/pins/orientation-symmetry/pins.py

MACEM_SUITE_ROOT overrides suite discovery for another checkout. No oracle files are written.
Calls the separate pinned py4DSTEM interpreter, and compiles read-only git-show app source for native keys/samples.
"""
import hashlib
import importlib.metadata as md
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from pymatgen.symmetry.groups import PointGroup
from scipy.spatial.transform import Rotation

KIT = Path(__file__).resolve().parents[3]
SUITE = Path(os.environ.get('MACEM_SUITE_ROOT', KIT.parent.parent if KIT.parent.name == '_worktrees' else KIT.parent))
OUT = KIT / 'Tests/DiffractoOrientationTests/Fixtures/symmetry-pins.json'
SYMBOLS = ['-1','2/m','mmm','4/m','4/mmm','-3','-3m','6/m','6/mmm','m-3','m-3m']
APP_PIN = 'c809679c'
SEED = 20261016

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def fold(d, ops):
    d=d/np.linalg.norm(d)
    best=d
    for op in ops:
        for sign in [-1,1]:
            v=sign*(op@d)
            delta=next((a-b for a,b in zip([-v[2],-v[0],-v[1]],[-best[2],-best[0],-best[1]]) if abs(a-b)>1e-13),0)
            if delta<0: best=v
    return best/np.linalg.norm(best)

def sample(count, ops):
    if count==0: return []
    selected=[]
    for d in np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=float):
        v=fold(d,ops)
        if not any(v@x>1-1e-13 for x in selected): selected.append(v)
    if count<=len(selected): return np.array(selected[:count]).tolist()
    n=max(192,count*16)
    golden=np.pi*(3-np.sqrt(5))
    candidates=[]
    for i in range(n):
        z=1-2*(i+.5)/n
        r=np.sqrt(max(0,1-z*z)); a=golden*i
        candidates.append(fold(np.array([r*np.cos(a),r*np.sin(a),z]),ops))
    candidates=np.array(candidates)
    distances=np.array([min(1-v@x for x in selected) for v in candidates])
    while len(selected)<count:
        best=int(np.argmax(distances))
        v=candidates[best];selected.append(v)
        distances=np.minimum(distances,1-candidates@v)
        distances[best]=-np.inf
    return np.array(selected).tolist()

def main():
    assert md.version('pymatgen') == '2026.9.24'
    assert md.version('pymatgen-core') == '2026.10.16'
    assert md.version('numpy') == '2.5.4'
    assert md.version('scipy') == '1.18.1'
    import pymatgen.symmetry.groups as groups
    clone = SUITE / 'References/DiffractoKit/pymatgen-core-ae9a834'
    core_pin='ae9a8348177b92fe33b214e54e61a3b8f4244b6e'
    assert subprocess.check_output(['git','-C',str(clone),'rev-parse','HEAD'],text=True).strip()==core_pin
    for f in ['groups.py','symm_data.json']:
        assert digest(Path(groups.__file__).with_name(f)) == digest(clone/'src/pymatgen/symmetry'/f)
    rng = np.random.default_rng(SEED)
    matrices = Rotation.random(32, random_state=rng).as_matrix()
    # Proper rotation matrices, rows lab axes in crystal coordinates after transpose convention.
    matrices = np.concatenate((matrices, [np.eye(3), Rotation.from_rotvec([0,np.pi/2,0]).as_matrix(),
        Rotation.from_rotvec([0,0,np.pi/4]).as_matrix(), Rotation.from_rotvec([np.pi,0,0]).as_matrix()]))
    halfturns = [Rotation.from_rotvec(angle*np.array(axis)/np.linalg.norm(axis)).as_matrix()
        for axis in [[np.cos(.1),0,np.sin(.1)],[np.cos(.5),0,np.sin(.5)],[1,2,3],[-2,3,1]]
        for angle in [np.pi,np.pi-1e-10,np.pi+1e-10]]
    matrices = np.concatenate((matrices,halfturns))
    directions = matrices[:,:,2]
    result = {'seed': SEED, 'python':sys.version.split()[0], 'versions': {k:md.version(k) for k in ['pymatgen','pymatgen-core','numpy','scipy']},
        'pymatgenCorePin':core_pin, 'pymatgenGroupsSHA256':digest(groups.__file__),
        'pymatgenSymmDataSHA256':digest(Path(groups.__file__).with_name('symm_data.json')), 'appPin':APP_PIN, 'groups':[]}
    for symbol in SYMBOLS:
        pg = PointGroup(symbol)
        frac = sorted([o.rotation_matrix.tolist() for o in pg.symmetry_ops])
        # pymatgen's tables use hexagonal fractional coordinates for trigonal/hexagonal groups.
        H = np.array([[1,0,0],[-.5,np.sqrt(3)/2,0],[0,0,1]]) if symbol in ['-3','-3m','6/m','6/mmm'] else np.eye(3)
        ops = [H.T @ np.array(m) @ np.linalg.inv(H.T) for m in frac if np.linalg.det(m)>0]
        cases=[]
        for M in matrices:
            candidates=[]
            for op in ops:
                C = op @ M
                q = Rotation.from_matrix(C.T).as_quat(scalar_first=True)
                if (next((x for x in q[1:] if abs(x)>1e-14),0)<0 if abs(q[0])<1e-14 else q[0]<0): q=-q
                candidates.append(([-q[0],*q[1:]],C,q[0]))
            best=candidates[0]
            for candidate in candidates[1:]:
                delta=next((a-b for a,b in zip(candidate[0],best[0]) if abs(a-b)>1e-13),0)
                if delta<0: best=candidate
            cases.append({'matrix':M.tolist(), 'reduced':best[1].tolist(), 'angle':float(2*np.arccos(np.clip(best[2],0,1)))})
        generic = symbol not in ['m-3m','6/mmm','2/m']
        result['groups'].append({'symbol':symbol, 'operators':[o.tolist() for o in ops], 'cases':cases,
            'foldedDirections':[fold(d,ops).tolist() for d in directions] if generic else [],
            'samples':[{'count':n,'directions':sample(n,ops)} for n in [0,1,2,3,8,16]] if generic else []})
    pyreq = {'symbols':SYMBOLS,'directions':directions.tolist(), 'matrices':matrices.tolist(),'groups':result['groups']}
    proc=subprocess.run([str(SUITE/'References/python/py4dstem/bin/python'),'-I','-B',str(KIT/'tools/pins/orientation-symmetry/pins-py4dstem.py')],
        input=json.dumps(pyreq),text=True,capture_output=True,check=True,
        env={**os.environ,'MPLCONFIGDIR':str(KIT/'.build/symmetry-mpl'),'XDG_CACHE_HOME':str(KIT/'.build/symmetry-cache')})
    result['py4dstem']=json.loads(proc.stdout)
    # Native keys and exact source samples. Only the extracted classes are compiled; the app is never built.
    source=subprocess.check_output(['git','-C',str(SUITE/'mac4DSTEM'),'--no-optional-locks','show',f'{APP_PIN}:mac4DSTEM/Core/Analysis/OrientationResult.swift'],text=True)
    result['appSourceSHA256']=hashlib.sha256(source.encode()).hexdigest()
    text='import Foundation\nimport simd\n'+source[source.index('/// Cubic crystal-symmetry'):source.index('/// Explicit symmetry policy')]
    text+=source[source.index('private extension simd_double3x3'):source.index('package nonisolated struct OrientationResult')]
    text=text.replace('package ','')
    literal=',\n'.join('SIMD3<Double>('+','.join(repr(float(x)) for x in v)+')' for v in directions)
    text+='\nlet directions = ['+literal+']\nvar records: [[String: Any]] = []\n'
    for symbol,cls in [('m-3m','CubicOrientationSymmetry'),('6/mmm','HexagonalOrientationSymmetry'),('2/m','MonoclinicOrientationSymmetry')]:
        text+='''
records.append(["symbol": "%s", "colors": directions.map { d in
 let c = %s.ipfColor(direction: d)
 return [c.x.bitPattern, c.y.bitPattern, c.z.bitPattern]
}, "directions": directions.map { d in
 let v = %s.reduceDirection(d); return [v.x,v.y,v.z]
}, "samples": [0,1,2,3,8,16].map { n in
 ["count":n, "directions":%s.sampleFundamentalZone(count:n).map { [$0.x,$0.y,$0.z] }] as [String:Any]
}])
'''%(symbol,cls,cls,cls)
    text+='let data = try JSONSerialization.data(withJSONObject: records, options: [.sortedKeys])\nprint(String(decoding:data,as:UTF8.self))\n'
    tempdir=KIT/'.build/symmetry-generator'
    tempdir.mkdir(parents=True,exist_ok=True)
    src=tempdir/'native.swift'; src.write_text(text)
    subprocess.run(['swiftc','-O','-module-cache-path',str(KIT/'.build/swift-module-cache'),str(src),'-o',str(tempdir/'native')],check=True)
    result['native']=json.loads(subprocess.check_output([str(tempdir/'native')]))
    result['nativeDirections']=directions.tolist()
    OUT.write_text(json.dumps(result, sort_keys=True, separators=(',',':'),allow_nan=False)+'\n')
    subprocess.run([sys.executable,'-I','-B',str(KIT/'tools/pins/orientation-symmetry/color-table.py')],check=True)
    print(f'Wrote {OUT}: {OUT.stat().st_size} bytes')

if __name__=='__main__': main()
