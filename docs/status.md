# Status

What is live, one table; then the hand-off. Rules are in `CLAUDE.md`; the registration in `docs/registration-slice-1.md`.

## Live (2026-10-15 early morning, end of the first session)

| Area | State | Evidence |
|---|---|---|
| Registration of slice 1 | Written before any code (`93fa321`); dated amendments in §8, the last after Gate B | `docs/registration-slice-1.md` |
| `DiffractoCore`: cells of the seven systems | `UnitCell`: lattices, metrics, d, g, the metric family; constructors per system | `UnitCellTests` (closed forms per system to 1e-12, the β″ anchors) |
| Symmetry operations and Laue groups | `SymmetryOperation` (the x,y,z grammar; singular rotations refused), `IntegerMatrix`, `PointSymmetry` (holohedries, Laue class, orbits, the canonical rule that reproduces the apps' names) | `PointSymmetryTests`, `SymmetryOperationTests` |
| Symmetry detection from sites | `SymmetryDetection`: the metric's holohedry (and a file's loop) tested on the sites; translations refined by least squares, exact from the file or on the k/24 grid only within their own uncertainty, re-verified; the set closed into a group | `SymmetryDetectionTests`, `RefuterRegressionTests` |
| CIF import | `CIFImport` / `CIFDocument`: the apps' parser and refusals, every system admitted, operations detected after a file's expansion, co-sited rows kept | `CIFImportTests` (β″ 22 atoms, 1448 / 419; the apps' fixtures; the refusals) |
| Reflections, absences, families | `ReflectionList`: per-axis tile with a size guard, the general reflection condition (exact for least-squares translations), macTEM's 1e-9 ring rule, Miller–Bravais names in hexagonal axes only | `CIFImportTests`, `PinsTests`, `RefuterRegressionTests` |
| Pins against pymatgen | 19 structures and 7 cells equal within the registered bounds (pymatgen-core 2026.9.23, spglib 2.8.0) | `PinsTests`; `tools/slice-1/pins.py` (6960 self-checks); `Fixtures/slice-1-pins.json` |
| Placeholder modules | `DiffractoScattering`, `DiffractoGeometry`, `DiffractoOrientation`: a scope statement each | their placeholder tests |
| Tests | 72 tests, 0 failures | `tools/test.sh` (release) at `16ce14d` |
| Gate | `tools/test.sh` green and stamped at `16ce14d` (72 tests, release); the frozen slice-1 commit `b4cfe17` was green with 56 | `.build/last-green` |
| Consumer check | `tools/consumer-check` builds in Debug with 0 warnings and passes 109 checks through the public API at `16ce14d` (skeleton: no negative controls yet) | `tools/consumer-check/run.sh` |
| Break-once proof | 29 of 29 mutants red on `16ce14d` (19 planned, 10 for the Gate B fixes) | `docs/slice-1-mutants-2026-10-14.md`, `tools/slice-1/mutants.py` |
| Gate B refuter | One refuter on `b4cfe17`: all 230 groups in standard settings agree; F1, F2, F4, F5 fixed; F3 and F6 recorded | `registration-slice-1.md` §8, 2026-10-15 |

Validation: none of this is validated on a real diffraction pattern; agreement with pymatgen and recovery of planted
cells are what is claimed (`registration-slice-1.md` §7).

## Hand-off

The next session starts from the suite board's `lanes/diffractokit` prompt. In order:

1. Gate B is closed: the fixes are in, gate and mutants green at `16ce14d`. A second refuter only on a named doubt
   (registration §6); none is named.
2. Slice 1b, registered before code: a lattice reduction (Niggli or Delaunay, then the conventional cell) so a P1 file
   in a non-conventional setting detects its full symmetry (Gate B F3), and a Cartesian position tolerance as an option
   (F6); faster detection for large symmetric P1 cells (F5, about N³ now).
3. Slice 2, registered before code: `DiffractoScattering` with py4DSTEM's Lobato table, F(hkl), |F|², the |F| > 1e-4 cut
   the apps' lists use; oracles pymatgen (X-ray F as a cross-check of the phase sums) and py4DSTEM's table; after the
   cut the families reproduce the apps' counts (Au 8, Si 10, Mg 24, β″ 417 / 1442).
4. The consumer-check's negative controls and `docs/consumers-0.1.0.md` come with the first tag.

What the apps will see when they bump (each app's own reviewed lane): a file's loop no longer decides the Laue class;
magnesium has a 25th symmetry-allowed ring before the |F| cut; the Si ring is named `{5 1 1} + {3 3 3}`; `x/2` and
singular operations are refused; a trigonal structure imports as -3 or -3m instead of being refused.

## Unverified

- No real CIF from a database (COD, ICSD) has been read yet; the fixtures are the apps', the β″ generator's and the
  refuter's. SKCIFParser (iRASPA) and Gemmi are named as parsing comparisons, not run.
- No performance measurement beyond the refuter's (a 40 Å cubic cell at 4 Å⁻¹: 17 million reflections, 6 s, 4.4 GB,
  now refused above 10 million).
