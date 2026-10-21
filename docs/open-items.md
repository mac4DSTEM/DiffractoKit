# Open items

A plain list, worst first. Each line: the item, the cause where it is known, what would close it.

- **`cif_import` is 27 times slower than pymatgen's `CifParser` on a 32-atom P1 supercell (about 67 ms against 2.5 ms; 351 ms
  against 2.5 ms in 0.9.0's `docs/benchmarks.md`), and `symmetry_detection` 5 times slower than spglib on the same cell (13 ms).** Not
  like for like: the kit's `CIFImport.read` also detects the symmetry; pymatgen's parser only parses. Cause, MEASURED on
  `tools/benchmarks/inputs/al_super2.cif` (1 536 operations, 0.9.0, in-process, median of runs): `read` 297 ms, of which
  `SymmetryDetection.detect` 287 ms, of which `closedSubset` alone 281 ms (the closure check, 95 %); the search for the operations
  and the parse are about 16 ms together. Not the conventional-cell step and not the search mode. The check composed every
  operation with every other (2.4 million pairs), each product built as a heap-allocated matrix and searched among the 32
  translations of its rotation. 0.9.1 (same operations, same tolerance semantics, `ClosureEquivalenceTests` against the old
  code): rotation ids and a product table, translations in a bucket list; `closedSubset` 58 ms, `read` 67 ms. What would refute
  the cause: a profile with `closedSubset` stubbed out would still show a large `read`; it does not (the rest is about 10 ms).
  Still open: 58 ms of closure is 2.4 million table and bucket lookups. Close with a generator-based closure (check the product
  of each operation with a generating set, not with every operation), which needs an argument that it drops the same operations
  on the pseudo-symmetric sets the split-atom tests plant; or skip the check when the file's loop already lists a closed set.
- **`reflections(kMax:)` is 2.3 to 6.6 times slower than `get_points_in_sphere`.** Not like for like: the kit also gives d, the
  absence flag, the rings and their orbits. Close by timing the kit's own enumeration without the ring grouping if it matters.
- **Phase matching.** Defaults unchanged. (a) `minimumMatrixEvidence` (new, off) keeps vacuum and stray-peak positions out of the
  matrix class; (b) `completenessAwareCrossPhaseRanking` stays off: it loses at small noise and wins near the pair radius
  (`PhaseDiagnosisTests`). Both need one labelled scan to be switched on: the owner's study scan with the phases marked.
- **ACOM.** 7 of 200 off-grid self-recoveries land 36 or 60 degrees away (zone axes whose centrosymmetric pattern has a symmetry
  of its own; not separable on spots alone). "Zone axis misses up to 12.8 degrees beyond the bank sampling" (mac4DSTEM's note)
  is not reproduced here: within 2.9 degrees of a template the nearest-bank gap is at most 1.6 degrees at 600 axes; a pattern
  between templates in a sparse corner of the zone is untested. Metal is not faster than all CPU cores up to 1 024 patterns
  at 600 axes; 16 384 patterns and 2 400 axes are the nightly's.
- **CIF symmetry.** Still unchecked: a group named by Hall symbol only, and a partial list whose count equals the primitive
  group's order but whose operations are the wrong ones (a count cannot see that).
- **Python pins.** `References/python/orix` and `py4dstem` hold pymatgen-core 2026.10.16 under a `pymatgen` 2026.9.24 metapackage
  (dist-info of both); `CLAUDE.md` and `docs/parity.md` say 2026.9.24, the older fixtures record pymatgen-core 2026.9.23 (a
  previous venv). `cif-symmetry-pins.json` records `pymatgen` 2026.9.24 (a literal in `pins.py`) and was generated on this
  install, so on pymatgen-core 2026.10.16: regenerating it on 2026-10-21 gave a byte-identical file. The other fixtures were not
  regenerated; re-generate and re-pin them together at the next fixture change.
