# tools/pins/cells: the pin generator

`pins.py` writes `Tests/DiffractoCoreTests/Fixtures/pymatgen-cells-pins.json` before any Swift exists. The pins are what pymatgen says; the
Swift tests read them and the planted truth. A difference is a failure of the port, never a tolerance to widen.

```
../References/python/pymatgen/bin/python -I -B tools/pins/cells/pins.py      # from any directory, about 5 s, deterministic
```

The oracle is pymatgen 2026.9.24 with pymatgen-core 2026.10.16, spglib 2.8.0, NumPy 2.5.4, Python 3.12 (read at run
time into `versions`). The beta'' CIF is written first by `tools/crystal-structures/make_beta_double_prime.py` into the
gitignored `References/crystal-structures/`; its sha256 and full text are in the pins (`cif_text`) so the Swift tests
rebuild it byte for byte. Every other CIF is a fixture in `Tests/DiffractoCoreTests/Fixtures/` whose first line is
`# provenance: ...`.

## What it writes

`cells`: the seven parameter cells (monoclinic first): a, b, c, angles, the lattice matrix (`vesta=True`:
a along x, b in xy), metric, volume, reciprocal rows, and |g|, d for every (h, k, l) in -3 ... 3 except 0.
`structures`: beta'' from the CIF, beta'' as a P1 list of its 22 expanded sites (symmetry detected, not read), and the
fixtures with planted truth (plus a few extras, named below). Per structure: cell, matrix, metric, volume, reciprocal, the expanded
sites (Z, fractional in [0, 1), occupancy), the declared IT number, spglib's number and point-group symbol, the Laue class,
every operation (integer rotation, translation, `x,y,z` string), the Laue rotations {R} and {-R}, every reflection at
k_max 1.2 Å⁻¹ (g, d, allowed, canonical form; the orbit only on the first 24 reflections of each structure, to keep the file small) and the families on the allowed reflections; `families_all` counts the
families on every reflection and gives the first 12 names; beta'' also has counts at 1.6 and 2.0.

## Definitions applied (written before the code, not from any Swift)

- Lattice: `Lattice.from_parameters(..., vesta=True)` rebuilt from the cell the CIF states; `g = h b1 + k b2 + l b3`
  with rows b from `reciprocal_lattice_crystallographic` (asserted a_i . b_j = delta_ij and B = G^-1 A); d = 1/|g|
  (asserted equal to `Lattice.d_hkl` to 1e-12 relative for every reflection).
- Operations: `SpacegroupAnalyzer(structure, symprec=1e-4).get_symmetry_operations(cartesian=False)`; rotations
  asserted integer; translations mod 1 in [0, 1), snapped to the nearest k/24 (asserted within 1e-9; pymatgen has
  already rounded spglib's to `Fraction.limit_denominator(1000)`). The set is asserted distinct, closed under
  composition, and a permutation of the sites.
- Laue class: from the rotation set {±R} by its order and its proper rotations of order 4 and 6, asserted equal to the
  32 -> 11 table applied to pymatgen's point-group symbol. "A 4-fold / a 6-fold" is read as a proper rotation (det +1): read
  over every element, -3 (order 6) makes -3m, m-3 and the -3m of CdI2 come out as 6/m or 6/mmm (the run prints
  the literal reading beside the right one).
- Orbit: `{hkl @ R : R in the Laue rotations}` (row vector times the matrix that acts on fractional coordinates, no
  transpose); canonical form = the member that is greatest under (number of non-negative indices, (h, k, l)).
  Asserted: the orbit is closed, the canonical form is constant on an orbit and differs between orbits.
- Allowed: absent iff some operation has `hkl @ R == hkl` and `hkl . t` is not an integer within 1e-9.
- Reflections: `get_points_in_sphere` on the crystallographic reciprocal lattice, asserted equal to the brute-force box
  |h| <= ceil(k_max |a_i|) and its distances equal to |g| to 1e-12. Order: increasing |g|; ties (within 1e-12 relative of the
  cluster's first) by descending (h, k, l).
- Families: sweep of the allowed reflections; a reflection joins while g <= g0 (1 + 1e-9); g = smallest member;
  members descending; orbits = distinct canonical forms descending; name = braces joined by ` + `, four Miller-Bravais indices
  (i = -(h + k) third) for the Laue classes -3, -3m, 6/m, 6/mmm. So the Si ring is `{5 1 1} + {3 3 3}` (the descending order
  macTEM's `Set.sorted(by: >)` gives), not the `{3 3 3} + {5 1 1}` of the first written prediction.
- Truth checks run in the same script (they stop it if they fail): the closed forms for all seven cells against the metric
  tensor, the macTEM anchors of the monoclinic cell, the reflection conditions and counts that are group theory.

## Notes the Swift side needs

- pymatgen's `CifParser` rounds any fractional coordinate within 1e-4 of 1/3 or 2/3 to the exact fraction before expanding
  (`io/cif.py` around line 610). The pinned `sites` therefore carry exactly 1/3 where a file writes 0.3333333 (magnesium,
  fluorapatite, WS2, CdI2); for those the pins add `sites_as_written` (parsed with `frac_tolerance=0`) and `site_snap_max`.
  The symmetry is the same either way.
- Translations: a detected translation from coordinates written to 6-7 digits carries 1e-7 noise; the absence test at 1e-9
  needs it snapped (to k/24) first, as the pins do.
- `file_operations` says whether the CIF's own loop equals the detected set (gold with 4 centring operations and the 4-operation
  tetragonal file do not: their loops hold fewer operations than the atoms' symmetry).
- Fixtures changed from the apps' test texts (each carries a `MODIFIED` comment): `calcite_R-3c_hex` (36 operations, the seed lists 3),
  `alpha_quartz_P3221` (operations 4-6 in the setting whose 2-fold axes pass z = 0, the seed's give 12 atoms), `pyrite_Pa-3` (24
  operations, the seed lists 12). The texts as seeded are `alpha_quartz_P3221_as_seeded` and `pyrite_Pa-3_as_seeded`.

## Tetragonal settings

```
../References/python/pymatgen/bin/python -I -B tools/pins/cells/tetragonal-pins.py
```

`tetragonal-pins.py`: pymatgen 2026.9.24 / spglib 2.8.0 (symprec 1e-4) on rutile, anatase and beta-Sn in three settings and two
controls; writes `Tests/DiffractoCoreTests/Fixtures/tetragonal-pins.json`. The kit equals spglib in every setting (2026-10-15).
