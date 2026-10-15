# DiffractoKit 0.1.1: what mac4DSTEM and macTEM change to re-pin

Written 2026-10-15 by DiffractoKit's session, for mac4DSTEM's and macTEM's sessions. 0.1.1 is slice 1b: the symmetry of a
structure in any cell setting, a Cartesian tolerance, and fast detection. The tag is commit 4876b18 (it also carries the CI
type-check fix). Registered in `docs/registration-slice-1b.md` (with its §6 amendments: the port, held back from 0.1.0, restored
for 0.1.1). The tag is the owner's; until it exists nothing here can be pinned. Read `docs/consumers-0.1.0.md` first: this note
lists only what changed since 0.1.0.

## In one table

| Who | Deletes or changes | Why |
|---|---|---|
| both | a file in a non-conventional setting now reads its full symmetry | `CrystalStructure.detectingSymmetry` detects in the conventional cell by default (`conventional: true`) |
| both | nothing may assume `operations` are in the conventional cell | they stay in the FILE's setting; the conventional structure is `conventionalSetting` |
| mac4DSTEM | the metric test in its Laue class to ACOM family mapping may reduce first | `UnitCell.conventionalCell()` and `niggliReduced()` (below); the mapping still tests `metricFamily()` |
| macTEM | a goldens-of-record on a monoclinic cell with a = c or a P1 file of a centred structure | the number moves (below); `conventional: false` keeps the file setting's reading |
| both | nothing is required for the new tolerance and search parameters | defaults are slice 1's (`.fractional(5e-3)`, cell-list search with the same operations bit for bit) |

## What numbers move

1. **A non-conventional setting now reads the symmetry of its lattice.** Before 1b the operations were detected in the cell as
   given, so a P1 file of a centred structure, a primitive cell (fcc or bcc in its 60° or rhombohedral primitive cell, a
   C-centred primitive cell) or a cell whose metric is higher than its declared system read only the symmetry its axes happened to
   show. `detectingSymmetry` (with `conventional: true`, the default) now carries the sites into the conventional cell, detects
   there by slice 1's rule, and reports back: the Laue class, the operation count, the absences and the merged rings are the
   full ones. The truths (registration §3, all met): a gold primitive cell is m-3m with 192 operations in the conventional cell
   and 8 rings, bcc iron m-3m with 96, β″ (C-centred primitive) 2/m with 8 operations and 417 rings, 1442 kept after the cut;
   rutile with its unique axis a is 4/mmm with 16 operations; magnesium with γ = 60° is 6/mmm with 24; corundum in
   rhombohedral axes is -3m, 36 operations in hexagonal axes, its ring names in Miller-Bravais form there.
2. **The file's own setting is reported, never relabelled.** `CrystalStructure.operations`, `laueGroup`, `laueClass`, the
   reflections and the rings stay in the FILE's cell (found through the conventional cell and conjugated back: R' = P R P^-1,
   t' = P t; reflections h_conv = h_file P). `CrystalStructure.setting` says which cell the numbers are in (`.file` for
   everything a caller builds or imports; `.conventional` for `conventionalSetting.structure`).
   `CrystalStructure.conventionalSetting: ConventionalSetting?` carries the same structure in the conventional cell (`structure`,
   `conventionalCell`, `cell`, `bravais`, `transformation`, `sites`, `operations`, `laueGroup`, `laueClass`,
   `translationUncertainty`) and maps between the two: `conventionalIndices(of:)`, `fileIndices(ofConventional:)`,
   `fileOperation(of:)`, `conventionalOperation(of:)`. It is nil when the operations were given (a CIF's own loop) or when
   `conventional: false`.
3. **A worked example of a number that moved.** The monoclinic cell with a = c and cos β = -1/4 is metrically C-centred
   orthorhombic: it now reads mmm, and {1 0 0} and {0 0 1} are one orbit. The kit's own `RingIndexingTests` preset (`mSq`) is the
   case; macTEM's accidental-coincidence reading is the file setting's 2/m, so the kit's test detects it with
   `conventional: false` (and also holds the merged orbit under full detection). No other existing test moved.
4. **Large symmetric supercells detect in about a second.** The permutation test uses a cell list (bins of the position
   tolerance, O(1) per image). A 5x5x5 gold supercell (24000 operations) takes 1.1 s (slice 1: about 30 s); a 500-atom
   symmetric P1 cell 0.008 s. The operations and translations equal slice 1's bit for bit (69 fixtures, against
   `SymmetryDetection.Search.bruteForce`, the CPU reference path kept).

Nothing else moves: the closed forms, scattering, the ring fit and the CIF reader are 0.1.0's bit for bit. A structure whose
operations are given (an imported CIF's loop) keeps its path; `conventionalSetting` is nil there.

## What is new (additive)

- `UnitCell.niggliReduced(relativeEpsilon:) -> ReducedCell` (`cell`, `transformation: IntegerMatrix` with det +-1, `epsilon`,
  `converged`, `steps`; `UnitCell.niggliRelativeEpsilon` = 1e-5; Krivy-Gruber with the Grosse-Kunstleve, Sauter and Adams
  variant).
- `UnitCell.conventionalCell(tolerance:) -> ConventionalCell` (`cell`, `bravais: BravaisLattice`, `transformation`, `reduced`,
  `centringTranslations`, `uniqueAxis`, `setting`, `centringMultiplicity`; `UnitCell.conventionalCellTolerance` = 1e-9).
  `BravaisLattice` is aP, mP, mC, oP, oC, oI, oF, tP, tI, hR, hP, cP, cI, cF (`centring`, `holohedryOrder`,
  `centringMultiplicity`). Setting rules the code settled: monoclinic unique axis b with C centring and beta >= 90 deg, hR in
  obverse hexagonal axes, and among equivalent bases the one closest to the identity; epsilon = 1e-5 V^(2/3). The 44 lattice
  characters are not tabulated: the cell is derived from the lattice's integer automorphisms on the Niggli basis (entries -1, 0,
  1), which is equivalent.
- `CrystalStructure.detectingSymmetry(..., tolerance:, conventional:, search:)`: `tolerance` is
  `SymmetryDetection.Tolerance`, `conventional` defaults to true, `search` is `SymmetryDetection.Search`.
- `SymmetryDetection.Tolerance`: `.fractional(5e-3)` (the default, slice 1's, `.default`) or `.cartesian(Å)` (spglib's
  sense, the minimum-image distance in Å); `.proposedCartesian` = `.cartesian(1e-3)` (10 x spglib's symprec 1e-4 Å). The
  default stays fractional until a measured comparison is registered; a pseudo-symmetry of 0.013 Å in a 50 Å cell reads as
  symmetry at the fractional default and not at `.cartesian(1e-3)`.
- `SymmetryDetection.Search`: `.cellList` (default) and `.bruteForce`.
- `CrystalStructure.Setting` (`.file`, `.conventional`), `.setting`, `.conventionalSetting`; `ConventionalSetting` as above.

## Evidence

Oracles: spglib 2.8.0 and pymatgen-core 2026.9.23 (`References/pymatgen-venv`). Niggli reduction against both on 2827 cells,
worst 3.5e-11 (the oracles agree to 2.8e-13); Bravais names 2800 of 2800; operation counts and Laue classes equal spglib's in both
settings on 27 structures; 4332 rings to 1e-8 in g. Every new test broken once: mutants 27 of 27 red.

**Gate B on 4876b18 (one refuter): no wrong number or class on realistic inputs.** 4500 cells at the boundaries between
lattice characters, 11,300 Bravais names, 932 and 1200 structures, 3255 CIFs, and the cell list bit equal to the scan on 4800
cases. Four findings, recorded here because the apps read the cases:
- **F1, a cell sheared by more than about 300 to 1000 units was silently unreduced** (a cell far outside any real file). Fixed
  after the tag, in 0.2.0, by multi-unit shears and a non-silent fallback.
- **F2, a conventional monoclinic cell may be relabelled a <-> c** by the shortest-pair rule. Both settings are reported, so the
  transformation says what happened; an app compares in the setting it names, not by axis letters.
- **F3, `metricAngleToleranceDeg` does not act on the conventional route.** A bcc primitive cell written with 109.4712 deg
  (five decimals short of the exact angle) reads hR unless `metricLengthTolerance` is at least 1e-5, where spglib's default reads
  Im-3m. A consumer reading rounded primitive cells passes `metricLengthTolerance: 1e-5` or larger.
- **F4, `.cartesian(t)` acts like about t/2 of spglib's symprec** (measured: operations were lost at a displacement of 0.5 to 0.9 t against 1.0 to 1.1 t for spglib, in 15 of 15 sweeps; the composition check uses the same bound). Compare
  tolerances by that factor; the default is not affected.

## The consumer check

`tools/consumer-check` (a separate package, the public API only, Debug, 0 warnings) passes 537 checks (404 at 0.1.0), each new
1b symbol with a real call and a checked result, and its negative controls (`run.sh --controls`, `controls.py`) cover 1b's
symbols: each is reverted from `public` to `package` in a copy and the consumer build must fail with the diagnostic naming it.

## mac4DSTEM

Re-pins; reads full symmetry for P1 and primitive files with no code change. Its Laue class to ACOM family mapping (0.1.0's
note: m-3m with `metricFamily() == .cubic`, 6/mmm with `.hexagonal`, 2/m with `.monoclinic(uniqueAxis: .b)`, else `.identity`)
must still test the metric family, because the file setting's frame is what the orientation maps use: an fcc structure in its
primitive 60-degree cell now reads m-3m with 192 conventional operations, and its Cartesian frame is still not the cubic axes.
With 1b a consumer can reduce first: `cell.conventionalCell()` gives the cell, `bravais` and `transformation` that make the frame
cubic, hexagonal or monoclinic, and `conventionalSetting.structure` is the structure in it; whether the ACOM lane maps in the
conventional frame (and converts the orientations back by the transformation) or keeps the file frame and the metric test is its
own decision, in its lane.

## macTEM

Re-pins; a CIF or cell in a non-conventional setting merges rings that were split before (the primitive cells, the accidental
a = c monoclinic of item 3), and its goldens of such cases move: each is a registered change, not a regression to wave through.
Re-measure the ring multiplicities of the primitive-cell cases. The accidental-coincidence case keeps its reading by calling
`detectingSymmetry(..., conventional: false)`. Index display: a ring named in the file setting stays in the file setting;
`conventionalSetting?.conventionalIndices(of:)` gives the conventional name where the room wants both.
