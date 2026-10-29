"""pymatgen standard conventional/primitive cells of committed MP-shaped JSON.
Run with References/python/pymatgen/bin/python -I -B under a memory cap.
The input records are copied from pinned app tests, NOT queried from MP in tests.
The Al v3 envelope is the app's recorded fetch; Si/Mg are synthetic app test inputs.
No expected value is calculated by the Swift port. Symprec 1e-4 A, angle tolerance 5 deg.
"""
from pathlib import Path
from importlib.metadata import version
import json, platform
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import spglib
assert version('pymatgen') == '2026.9.24'
assert version('pymatgen-core') == '2026.9.23'
assert np.__version__ == '2.5.4'
assert spglib.__version__ == '2.8.0'
root=Path(__file__).resolve().parents[3]/'Tests/DiffractoCoreTests/Fixtures/mp-import'

def result(s):
    return dict(lattice=s.lattice.matrix.tolist(),sites=[dict(species={str(k):v for k,v in site.species.items()},abc=site.frac_coords.tolist()) for site in s])

cases=[]
for path in sorted(root.glob('input-*.json')):
    raw=json.loads(path.read_text()); doc=raw['data'][0]; st=doc['structure']; l=st['lattice']
    lattice=Lattice(l['matrix']) if 'matrix' in l else Lattice.from_parameters(l['a'],l['b'],l['c'],l['alpha'],l['beta'],l['gamma'],vesta=True)
    s=Structure(lattice,[{v['element']:v['occu'] for v in p['species']} for p in st['sites']],[p['abc'] for p in st['sites']])
    analyzer=SpacegroupAnalyzer(s,symprec=1e-4,angle_tolerance=5)
    cases.append(dict(input=path.name,spaceGroupNumber=analyzer.get_space_group_number(),spaceGroupSymbol=analyzer.get_space_group_symbol(),
                      conventional=result(analyzer.get_conventional_standard_structure()),primitive=result(analyzer.get_primitive_standard_structure())))
out=dict(python=platform.python_version(),pymatgen=version('pymatgen'),pymatgenCore=version('pymatgen-core'),spglib=spglib.__version__,numpy=np.__version__,symprec=1e-4,angleTolerance=5,cases=cases)
(root/'pymatgen-pins.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
print(f'{len(cases)} structures pinned')
for c in cases:
 print(c['input'],c['spaceGroupSymbol'],c['spaceGroupNumber'],len(c['conventional']['sites']),len(c['primitive']['sites']))
