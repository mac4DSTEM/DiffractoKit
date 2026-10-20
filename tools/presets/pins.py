"""pymatgen / spglib pins for the preset crystals (DiffractoKit 0.6.0, CrystalPresetsTests).

  ../tools/memcap.sh -- ../References/DiffractoKit/pymatgen-venv/bin/python -I -B tools/presets/pins.py

Every preset is built here independently of the Swift code, from the published lattice constants and Wyckoff positions with
`Structure.from_spacegroup` (so the sites are spglib/pymatgen's own orbit expansion, not a copy of the hand lists):
the cubic and hexagonal elements from their space groups, beta'' from C2/m (No. 12, unique axis b) and 2H-WS2 from
P6_3/mmc (No. 194). Recorded per preset: the expanded sites (fractional, wrapped to [0,1)), the cell volume, the
space-group number and symbol from SpacegroupAnalyzer (symprec 1e-5), the number of space-group operations spglib finds,
the formula, and d(hkl) for a list of reflections (pymatgen Lattice.d_hkl). Writes
Tests/DiffractoCoreTests/Fixtures/presets-pins.json.
"""
from importlib.metadata import version
import json
from pathlib import Path
import platform
import numpy as np
import spglib
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

assert version('pymatgen') == '2026.9.24' and version('pymatgen-core') == '2026.9.23'
assert np.__version__ == '2.5.4' and spglib.__version__ == '2.8.0'
OUT = Path(__file__).resolve().parents[2] / 'Tests/DiffractoCoreTests/Fixtures/presets-pins.json'
HKL = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (0, 0, 2), (1, 0, 1), (1, 0, 2), (2, 1, 1), (0, 2, 0), (-2, 0, 1), (1, -1, 3), (2, 2, 2), (0, 0, 4)]

def cubic(sg, a, sym, coords):
    return Structure.from_spacegroup(sg, Lattice.cubic(a), [sym], coords)

PRESETS = {
    'au_fcc': cubic('Fm-3m', 4.0782, 'Au', [[0, 0, 0]]),
    'al_fcc': cubic('Fm-3m', 4.0495, 'Al', [[0, 0, 0]]),
    'ni_fcc': cubic('Fm-3m', 3.5240, 'Ni', [[0, 0, 0]]),
    'cu_fcc': cubic('Fm-3m', 3.6149, 'Cu', [[0, 0, 0]]),
    'fe_bcc': cubic('Im-3m', 2.8665, 'Fe', [[0, 0, 0]]),
    'si_diamond': cubic('Fd-3m:1', 5.4309, 'Si', [[0, 0, 0]]),
    'sc_cubic': cubic('Pm-3m', 3.0, 'Po', [[0, 0, 0]]),
    # 2c at (2/3, 1/3, 1/4), then the origin moved to the first atom: the applications' (0,0,0), (2/3,1/3,1/2) setting.
    'mg_hcp': Structure.from_spacegroup('P6_3/mmc', Lattice.hexagonal(3.2094, 5.2108), ['Mg'], [[2 / 3, 1 / 3, 0.25]]),
    'ws2_2h': Structure.from_spacegroup('P6_3/mmc', Lattice.hexagonal(3.1532, 12.323), ['W', 'S'],
                                        [[1 / 3, 2 / 3, 0.25], [1 / 3, 2 / 3, 0.6225]]),
    'beta_double_prime': Structure.from_spacegroup(
        'C2/m', Lattice.from_parameters(15.16, 4.05, 6.74, 90, 105.3, 90), ['Mg', 'Mg', 'Mg', 'Si', 'Si', 'Si'],
        [[0, 0, 0], [0.3459, 0, 0.089], [0.430, 0, 0.652], [0.0565, 0, 0.649], [0.1885, 0, 0.224], [0.2171, 0, 0.617]]),
}

def main():
    PRESETS['mg_hcp'].translate_sites(list(range(2)), [-2 / 3, -1 / 3, -0.25], frac_coords=True, to_unit_cell=True)
    out = {'python': platform.python_version(), 'pymatgen': version('pymatgen'), 'pymatgenCore': version('pymatgen-core'),
           'spglib': spglib.__version__, 'numpy': np.__version__, 'symprec': 1e-5, 'hkl': HKL, 'presets': {}}
    for key, s in PRESETS.items():
        sga = SpacegroupAnalyzer(s, symprec=1e-5)
        ops = sga.get_symmetry_operations()
        sites = [{'z': site.specie.Z, 'abc': [float(x % 1.0) for x in site.frac_coords]} for site in s]
        out['presets'][key] = {
            'lattice': {'a': s.lattice.a, 'b': s.lattice.b, 'c': s.lattice.c, 'alpha': s.lattice.alpha,
                        'beta': s.lattice.beta, 'gamma': s.lattice.gamma},
            'volume': s.volume, 'spaceGroupNumber': sga.get_space_group_number(), 'spaceGroupSymbol': sga.get_space_group_symbol(),
            'pointGroup': sga.get_point_group_symbol(), 'crystalSystem': sga.get_crystal_system(),
            'operationCount': len(ops), 'siteCount': len(sites), 'formula': s.composition.formula,
            'sites': sites, 'd': [float(s.lattice.d_hkl(h)) for h in HKL]}
        print(key, out['presets'][key]['spaceGroupSymbol'], out['presets'][key]['spaceGroupNumber'], len(sites), len(ops))
    OUT.write_text(json.dumps(out, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')

if __name__ == '__main__':
    main()
