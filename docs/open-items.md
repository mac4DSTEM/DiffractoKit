# Open items

A plain list, worst first. Each line: the item, the cause where it is known, what would close it.

- **`symmetry_detection` is 28 times slower than pymatgen's spglib route on a 32-atom P1 supercell (361 ms against 13 ms;
  the same on a 4-atom cell, `docs/benchmarks.md`).** Cause (read from the code, not profiled): `SymmetryDetection.closedSubset`
  composes every operation with every other (1 536 operations: 2.4 million compositions, each a linear search of the 32 translations
  of its rotation, about 75 million comparisons). Not the conventional-cell step and not the search mode (all three timed at 300 ms).
  Close with a hashed closure check or a generator-based one, keeping the same operations (the fixtures in `CIFImportTests` pin them).
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
- **Python pins.** `refs/orix` holds pymatgen 2026.10.16; `docs/parity.md` pins 2026.9.24. Timing only here; re-pin before the
  next fixture is generated from it.
