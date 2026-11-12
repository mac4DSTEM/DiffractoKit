# Status

DiffractoKit 0.12.0 (2026-11-06), on ImagoKit `"0.2.0"..<"1.0.0"`.

Live:

- Cells in the seven systems, Niggli reduction, the conventional cell.
- CIF import with symmetry detection and a completeness check by IT number, by symbol or by Hall symbol (a Hall-only CIF
  is expanded with its group; an operation list that is not its space group is refused, and one that disagrees with the
  Hall symbol warns); space-group-free symmetry from the atoms; reflections and rings.
- Kinematic structure factors (py4DSTEM's scattering-factor table, bit for bit).
- Pattern geometry: ring and ellipse fits, calibration, zone axes, indexing; the ring-sequence assignment of a known crystal,
  Δq from several rings, and the single-pattern beam centre (`PatternCentre.find` refuses a blocked beam; `FriedelCentre`
  finds the centre of a pattern with one), radial profile and spot refinement.
- The ACOM matcher (CPU reference, Metal option, optional in-plane refinement) and vector-matched phase maps.
- Orientation symmetry and IPF keys, misorientation between two groups, `Orientation` from Euler angles with `dot`,
  `dotOuter` and `distanceMatrix`, orientation-relationship variants (beta'' in Al: 12), orientation-map planes.

Since 0.10.0: 0.11.0 (2026-10-30) dropped the unused public API; 0.12.0 (2026-11-01) added the items above, made the ACOM
Metal path 1.9 to 3.2 times faster and the ring profile fit stop runs that cannot improve, and retired the apps-parity dumps.
Details and the migration notes: `CHANGELOG.md`.

## In progress, not wired: phase orientations from the structure alone (#81)

The aim: a crystal structure in, and the scan's phases and their orientations out, with no zone axis typed. Five steps are
merged as kit functions; none is used by `PhaseVectorMatcher.map`. State of 2026-11-11.

| Step | PR | What it is |
|---|---|---|
| 1 | #82 | The matrix orientation of every position, its spots predicted and removed (`matrixOrientations`, `removingMatrixSpots`) |
| 2 | #84 | One orientation per grain (`regularizingMatrixOrientations`) and the guarded free search of the residual (`searchOrientations`) |
| 3 | #85 | Relationships learned from the strong positions, refined, and every position scored against their variants (`learnOrientationRelationships`, `classifyAgainstVariants`); grain overlaps in the removal |
| 4 | #86 | The two passes as one call, scored on the whole pattern, firm rule by default, neighbour rule as a switch (`classifyInTwoPasses`) |
| 5 | #88 | The rational relationships near a learned one, listed and never chosen (`rationalRelationships`) |
| guard | #95 | A displaced matrix pattern is never a phase (inside `classifyInTwoPasses`) |
| input | #96 | The detector's distortion measured on the matrix's spots, for a scan with no ring (`matrixSpotDistortion`) |
| one call | #99 | Steps 1 to 4 chained with the values they share, returning a `PhaseMap` (`mapFromStructures`); the call an application makes |

Where it stands against the bar (typed zone axes through `map`):

- Thronsen et al. dataset A, stride 3: 28 771 of 29 241 positions (98.39 %) with 103 matrix positions called a phase,
  against 28 858 (98.69 %) and 104. With the neighbour rule 98.51 % and 112. The positions that decide it show one vector of
  their own: the precipitate edges.
- The synthetic three-grain demo: relationship 0.03 degrees from the planted one, all nine precipitates found, nothing
  labelled outside them.
- Planted pairs (`TwoPassOtherMaterialsTests`):

| Pair | Relationship learned | Positions |
|---|---|---|
| alpha-Ti (hcp) + beta-Ti, Burgers | 2.65 degrees off, 12 variants for 6 | all right |
| alpha-Fe (bcc) + cementite, Bagaryatsky | 1.44 degrees off | all right |
| Al + beta'' (monoclinic), Andersen | exact | all right |
| Ni + gamma' Ni3Al, cube on cube | nothing learned | no gamma' found (known gap, pinned by the test) |

It is wired only when Thronsen A reaches the typed result on both numbers and nothing is worse on the demo or the pairs.

Next step: lattice-correspondence proposals in `searchOrientations`. For each position, besides the template proposals,
propose the precipitate orientations that put low-index precipitate lattice vectors onto matrix lattice vectors at that
position's matrix orientation (the free search does not take the matrix orientation today), and score them on the residual.
A phase on the matrix lattice (gamma') leaves only its superlattice spots after the removal and the template match never
proposes its orientation. The case to turn round is Ni + gamma' in `TwoPassOtherMaterialsTests`. Then before and after on
Thronsen A, the demo and the four pairs.

Parked: generating `map`'s entries (zone axis per variant, in-plane pair, slab) from a relationship. The typed Thronsen
entries carry hand-set values (the T1 zone [0 -4 1], a slab of 0.3 for theta'), and reproducing them by rule risks tuning
to that dataset.

Open beside it: the tolerances in Å⁻¹ should follow the disk size and the calibration (#87). Measurements, per step:
`docs/parity.md`.

## Tests

`DiffractoKitParityTests` holds the pinned comparisons against the Python references; the other test targets hold the
fast unit tests; `DiffractoKitTestSupport` holds the helpers both use.

What is checked against what: `docs/parity.md`. Timings against Python: `docs/benchmarks.md`. What is open: the GitHub
issues; known limits: `docs/parity.md`. Where it is going: `ROADMAP.md`.
