# DiffractoKit slice 1b: the symmetry of any cell setting, a Cartesian tolerance, fast detection (registration)

Registered 2026-10-15, about 03:00 CEST, before any code, by the DiffractoKit session (the owner's instruction of 00:55,
item 5: "if time remains: register slice 1b, registration only"). Tier 1: the Laue class, the absences and the ring names
are numbers and names a scientist reads. It closes three limits slice 1's Gate B recorded (registration-slice-1.md §7 and
§8, 2026-10-15: F3, F6 and the N³ import). No code is written in the session that registers it.

## 1. What slice 1b makes (DiffractoCore)

- **Lattice reduction.** `UnitCell.niggliReduced()` (Křivý & Gruber 1976, with the numerically stable variant of
  Grosse-Kunstleve, Sauter & Adams 2004, Acta Cryst. A60, 1–6, at a stated ε) returns the reduced cell and the integer
  change-of-basis matrix P (det ±1). `UnitCell.conventionalCell()` derives the conventional (centred) cell and its
  transformation from the reduced metric by the 44 lattice characters (International Tables A, §9.2, the Bravais-lattice
  table), at a stated tolerance; the result names the Bravais lattice (aP, mP, mC, oP, oC, oI, oF, tP, tI, hR, hP, cP, cI,
  cF) and the setting.
- **Symmetry of any setting.** `CrystalStructure.detectingSymmetry` gains `conventional: Bool = true`: the sites are
  transformed into the conventional cell (fractional x' = P⁻¹ x, wrapped, duplicates of the centring merged within the
  position tolerance), the operations detected there (slice 1's rule, unchanged), and the result reported in BOTH settings:
  the conventional structure (its operations, Laue class, absences and ring names) and the transformation back to the
  file's cell (operations conjugated: R' = P R P⁻¹, t' = P t; reflections h' = h P). A file's own setting is never
  relabelled silently: `CrystalStructure.setting` says which cell the numbers are in.
- **A Cartesian tolerance.** `SymmetryDetection.Tolerance` = `.fractional(5e-3)` (slice 1's, the default) or
  `.cartesian(Å)` (spglib's sense, default 1e-4 Å · 10 = 1e-3 Å proposed; the default stays fractional until a measured
  comparison is registered), applied to the minimum-image Cartesian distance of an image and its target.
- **Fast detection.** The permutation test by a cell list (bins of the position tolerance, or the Cartesian one), so each
  image looks up its candidate targets in O(1) instead of O(N): a 500-atom symmetric P1 cell in under 2 s (slice 1: 30 s),
  results identical to slice 1's (the same operations, the same translations bit for bit).

## 2. References and pins

- **spglib 2.8.0** (BSD-3-Clause, behind pymatgen-core 2026.9.23 in `References/pymatgen-venv`): `spglib.niggli_reduce`,
  `spglib.get_symmetry_dataset` (its `transformation_matrix`, `origin_shift`, `std_lattice`, `hall_number`) as the oracle
  for the reduced cell, the conventional cell and the operation set in the conventional setting. Agreement is agreement:
  spglib's own setting choices (origin choice, unique axis) are compared up to the registered equivalences below, and a
  difference in choice is not a defect.
- **pymatgen** `Lattice.get_niggli_reduced_lattice` (an independent implementation of the same algorithm) as a second
  oracle for the reduced cell.
- **International Tables A** (the 44 lattice characters, the Bravais lattices' reduced forms) as the definition; no table
  is copied, the conditions are written from the text.

Pins (`tools/slice-1b/pins.py`, written before the Swift): for every structure of slice 1's fixtures and the Gate B
refuter's non-conventional cells (fcc and bcc primitive cells, β″'s C-centred primitive cell, rutile with unique axis a,
magnesium with γ = 60°, corundum in rhombohedral axes, the 7× supercell, a triclinic cell whose reduced form is hexagonal),
plus 200 random cells of each Bravais lattice under a random unimodular transformation (|entries| ≤ 2): the reduced
cell's parameters (1e-10 relative), P (exact, up to the reduced cell's own symmetry: equal metric is the criterion), the
Bravais lattice name (exact), the conventional operation count and Laue class (exact), and the ring list after slice 2's
cut in the conventional setting (names, g 1e-12, multiplicities).

## 3. Truth

- A primitive cell of fcc gold, of bcc iron and of β″ (C-centred) detects m-3m (192 operations in the conventional cell),
  m-3m (96) and 2/m (8): the counts slice 1 gives for the conventional cells; the rings after the cut are the conventional
  cells' rings exactly (Au 8, β″ 417 / 1442).
- Rutile with unique axis a: 4/mmm, 16 operations; magnesium with γ = 60°: 6/mmm, 24; corundum in rhombohedral axes: -3m,
  in hexagonal axes for the conventional cell (36 operations), its ring names in Miller–Bravais form there.
- Every Bravais lattice's reduced form satisfies the Niggli conditions (the tests check the conditions, not only the
  oracle), and P is unimodular.
- The fast detection returns slice 1's operation set bit for bit on every slice-1 and slice-2 fixture.

## 4. Gates

Pins within their bounds; truth exact; every new test broken once (`tools/slice-1b/mutants.py`, slice 1's pattern);
one Gate B refuter on the frozen commit, asked to break the reduction near its ε (cells within 1e-6 of a boundary between
lattice characters) and the conventional transformation of centred cells with atoms on special positions; consumer-check
extended; `docs/consumers-0.0.4.md` says what moves for the apps (a P1 file in a non-conventional setting now gets its
full symmetry; a file's own setting is reported, not relabelled).

## 5. Not claimed

Settings that differ only by an origin choice (Fd-3m origin 1 or 2) are equivalent: the kit reports the one its reduction
reaches, with the shift. The Cartesian tolerance's default is not changed by this slice. No real CIF from a database is
read here (the parsing comparisons SKCIFParser and Gemmi stay a separate tier-2 item).

## 6. Amendments

### 2026-10-15, about 14:00 CEST, after the port (DiffractoKit's third session, lane E): what the code settled

- **The 44 lattice characters** are not tabulated in the code: the conventional cell is derived from the lattice's integer
  automorphisms on the Niggli basis (entries −1, 0, 1), which is equivalent (the characters are the classes those
  automorphisms define). Setting rules fixed by the port: monoclinic unique axis b with C centring and β ≥ 90°, hR in
  obverse hexagonal axes, and among equivalent bases the one closest to the identity. ε = 1e-5 · V^(2/3).
- **Both settings**: `CrystalStructure.operations`, the Laue group and the rings stay in the FILE's setting (found through
  the conventional cell, conjugated back); the conventional structure is `conventionalSetting.structure` with
  `.setting == .conventional`. The file's setting is never relabelled.
- **The rings against spglib** are compared at slice 1's level (before the |F| cut); the cut rings are compared in
  `DiffractoScatteringTests` (`ConventionalRingsTests`). Ring names are compared where the conventional cell has spglib's
  parameters (1294 compared, 872 skipped: β″, tetragonal a = c, the 7× supercell).
- **A number moves for one existing test** (the registration's "a non-conventional cell now gets its full symmetry"):
  `RingIndexingTests`' `mSq` preset, a monoclinic cell with a = c and cos β = −1/4, is metrically C-centred orthorhombic;
  with `conventional: true` it reads mmm and {1 0 0} and {0 0 1} are one orbit. macTEM's accidental-coincidence case is
  the file setting's 2/m, so that preset is detected with `conventional: false`, and the test now also holds the merged
  orbit under full detection. No other existing test moved.
- **Results**: Niggli reduction against spglib and pymatgen on 2827 cells, worst 3.5e-11 (the oracles agree to 2.8e-13);
  Bravais names 2800 / 2800; operation counts and Laue classes equal spglib's in both settings on 27 structures; 4332 rings
  to 1e-8 in g; the truths of §3 hold (Au primitive 192 operations and 8 rings, β″ primitive 417 rings and 1442 kept
  after the cut, rutile along a 4/mmm 16, Mg γ = 60° 6/mmm 24, corundum −3m 36 in hexagonal axes). Fast detection: a
  500-atom symmetric P1 cell in 0.008 s, a 5×5×5 gold supercell (24000 operations) in 1.1 s (slice 1: about 30 s), equal
  to the brute-force scan bit for bit on 69 fixtures. Mutants 27 / 27 red. The Cartesian default stays fractional (5e-3);
  `.cartesian(1e-3 Å)` is offered as `proposedCartesian`. What moves for the apps is in `docs/consumers-0.1.0.md` (slice 1b
  ships at 0.1.0, not 0.0.4).

### 2026-10-15, about 14:30 CEST: held back from 0.1.0

The owner decided (13:16) that the frozen 0.0.4 candidate is DiffractoKit 0.1.0 and slice 1b is 0.1.x. Its code
(commit bc99fa0, the port above) had already landed on `main` and has no Gate B refuter yet, so it is reverted on `main`
before the 0.1.0 tag and this registration and its amendment stay. The 0.1.x session restores it with
`git revert <the revert commit>` (or `git cherry-pick bc99fa0`), runs its one Gate B refuter (§4: the reduction near its
ε, centred cells with atoms on special positions), extends the consumer check (a prepared patch is named in the hand-off),
and writes what moves for the apps in `docs/consumers-0.1.1.md`.
