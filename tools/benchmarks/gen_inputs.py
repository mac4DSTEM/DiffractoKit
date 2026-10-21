"""Shared inputs of the DiffractoKit benchmark (both sides read these files).

Run: ../../../refs/orix/bin/python -I -B gen_inputs.py   (pymatgen, numpy)
Writes inputs/: cells.json (random cells), al_super{1,2,3}.cif (P1 fcc Al supercells), c2m.cif (C2/m with its operations),
ring_{128,256,512}.f32 (a clean elliptical Gaussian ring, row-major float32), rotations.f64 (random proper rotations, 3x3 row-major).
"""
import json, pathlib
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.cif import CifWriter

here = pathlib.Path(__file__).resolve().parent
out = here / "inputs"
out.mkdir(exist_ok=True)
rng = np.random.default_rng(20261021)

cells = []
for _ in range(10000):
    a, b, c = rng.uniform(3, 12, 3)
    al, be, ga = rng.uniform(60, 120, 3)
    cells.append([a, b, c, al, be, ga])
json.dump(cells, open(out / "cells.json", "w"))

al = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.0495), ["Al"], [[0, 0, 0]])
for n in (1, 2, 3):
    s = al.copy(); s.make_supercell([n, n, n])
    open(out / f"al_super{n}.cif", "w").write(str(CifWriter(s)))   # P1: every atom listed, no symmetry
c2m = Structure.from_spacegroup("C2/m", Lattice.monoclinic(11.6, 4.1, 6.5, 105.0), ["Mg", "Si", "Si"],
                                [[0.0, 0.0, 0.0], [0.19, 0.0, 0.22], [0.22, 0.0, 0.62]])
open(out / "c2m.cif", "w").write(str(CifWriter(c2m, symprec=1e-4)))

for n in (128, 256, 512):
    cx = cy = (n - 1) / 2 + 0.3
    a, b, th = 0.32 * n, 0.30 * n, 0.4
    y, x = np.mgrid[0:n, 0:n].astype(np.float64)
    dx, dy = x - cx, y - cy
    c, s = np.cos(th), np.sin(th)
    u = dx * c + dy * s; v = -dx * s + dy * c
    rho = np.sqrt((u / a) ** 2 + (v / b) ** 2)
    img = np.exp(-((rho - 1) * (a + b) / 2) ** 2 / (2 * 2.0 ** 2)).astype(np.float32)
    img.tofile(out / f"ring_{n}.f32")

q = rng.normal(size=(100000, 4)); q /= np.linalg.norm(q, axis=1)[:, None]
w, x, y, z = q.T
R = np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
              2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
              2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], axis=1)
R.astype(np.float64).tofile(out / "rotations.f64")
print("inputs written", sorted(p.name for p in out.iterdir()))
