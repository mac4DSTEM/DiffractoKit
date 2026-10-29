# tools/slice-2: the pin generators

Written before any Swift of slice 2 (DiffractoScattering) existed; the Swift tests
read the files these write. A difference is a failure of the port, never a tolerance to widen.

```
tools/slice-2/fetch-table.sh                                         # py4DSTEM f050d207 scattering_factors.txt into References/, sha256-locked
../References/python/orix/bin/python -I -B tools/slice-2/pins.py      # about 2 s, deterministic (no clock): writes slice-2-pins.json
```

- `pins.py` writes `Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json` (about 1.7 MB): the versions, the table's
  sha256 and the 103 x 10 table as `np.loadtxt` reads it, `f_e` at g2 in {0, 0.25, 1, 4} for Z = 1..103, 20 structures (the
  apps' presets Au, Al, Si, Mg, WS2; Cu3Au-occ and Weak-cut from macTEM a58c84c; ideal hcp; the extra non-centrosymmetric
  zincblende GaAs; beta'' and ten slice-1 CIF fixtures), each with every reflection above 1e-6 in |F|, the kept count, the
  rings after the cut, the ring counts on the operation-allowed and on all reflections, the operation-absent check; the
  mutation-moves; pymatgen's X-ray pattern per peak with the explicit half-space members; the
  `cut_edge` cell whose every |F| is exactly 1e-4; the registered predictions with the measured values.
- The oracle is NumPy written from the definitions (not from any Swift), the Laue group comes from pymatgen
  `SpacegroupAnalyzer` (symprec 1e-4), the ring rule and canonical forms are slice 1's (copied, not imported), and
  pymatgen's `XRDCalculator` is the independent phase sum. CIF structures are pinned on the sites AS WRITTEN
  (`frac_tolerance = 0`, what `CIFImport` keeps); `sites_snapped` and `snap_effect` carry pymatgen's default.
- Exit 0: every self-check and every registered predicted count holds. 1: a self-check failed. 3: the file was written but a
  predicted count differs (printed twice; the number is never edited to fit).
