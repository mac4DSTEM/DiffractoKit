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
  (`PhaseDiagnosisTests`). Both need one labelled scan to be switched on: a study scan with the phases marked.
- **ACOM.** 7 of 200 off-grid self-recoveries land 36 or 60 degrees away. Cause, measured: each is a coincidence twin of the
  truth (60 degrees about a <111> near <112> beams, 36.87 degrees about a <100> near <013> beams; the axis lies in the zone
  plane), whose kinematic pattern overlaps the truth's 0.95 to 1.0000; only the Ewald curvature term separates them, and
  py4DSTEM picks the twin in 7 of the 14 truth and twin patterns (`theLargeFailuresAreCoincidenceTwinsOfTheTruth`). Inherent
  to spot patterns at 200 kV, not a matcher defect; a consumer can flag it by scoring the twin of the winner. "Zone axis misses up to 12.8 degrees beyond the bank sampling" (mac4DSTEM's note)
  is not reproduced here: within 2.9 degrees of a template the nearest-bank gap is at most 1.6 degrees at 600 axes; a pattern
  between templates in a sparse corner of the zone is untested. Metal is not faster than all CPU cores up to 1 024 patterns
  at 600 axes; 16 384 patterns and 2 400 axes are not timed yet.
- **CIF symmetry.** A Hall symbol is read and a listed operation list is checked for its count, the cell metric, closure and
  the IT number's crystal class. Still unchecked: a closed list of the right class and count with the wrong translations
  (P4/m listed for P4_2/m, or another origin) when the file states no Hall symbol, since the kit has no table from a
  number or Hermann-Mauguin symbol to operations; and a full loop that disagrees with the Hall symbol beside it (the loop
  is used, as pymatgen does).
- **Beam centre for a blocked beam.** `PatternCentre.find` refuses a pattern whose beam is blocked; Friedel-pair symmetry would
  find it.
- **Ring fit.** Speed: `fitRing` at R 600 px with a band of ±90 px on 2048² took 24 s (py4DSTEM 13 s on one thread) because the
  profile optimiser ran to its 80-pass cap on a start stuck at the ring-width floor and after converging; with the two stops of the
  ring-fit speed fix it takes about 5 s, the same radius. Left: about 145 ms a pass at 680 000 pixels (the finite-difference Jacobian
  and the normal matrix: Accelerate, or a coarse-then-fine annulus). Neighbouring ring and small rings: one cause, the tail of the
  neighbouring ring inside the annulus; the bias follows the gap to the neighbour in ring sigmas (Au {111} with sigma 1.2 px: 0.03 to
  0.08 px at R 17 to 34 px, nothing from R 68), a single ring is exact at every radius, and py4DSTEM's `fit_ellipse_amorphous_ring`
  has the same bias in size (`RingFitNeighbourTailTests`). Shared by the reference, so not changed. Close with a profile model that
  has a term for the neighbouring ring (a new method, not a port), or keep the annulus inside the gap.
- **orix API gaps.** The orientation (Rodrigues) fundamental zone, and `dot` / `dot_outer` / `get_distance_matrix` across two
  different groups (`from_euler`, `dot`, `dot_outer` and `get_distance_matrix` for one group are held). orix's
  `Misorientation.reduce()` misses the smallest angle when the two groups' product set is not a group (432 with 622, 23 with
  3): it builds its zone from start . end while its equivalents need end . start (cause and pins in `Misorientation.swift`);
  the kit returns the smallest angle. Not reported upstream yet.
- **Not yet used by a known consumer** (0.10.0): `PatternCentre.find`, `PatternRadialProfile.about`, `LinearScale.fit` and
  `Misorientation`. Revisit at 0.12.
