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

Tests: `DiffractoKitParityTests` holds the pinned comparisons against the Python references; the other test targets hold the
fast unit tests; `DiffractoKitTestSupport` holds the helpers both use.

What is checked against what: `docs/parity.md`. Timings against Python: `docs/benchmarks.md`. What is open: the GitHub
issues; known limits: `docs/parity.md`. Where it is going: `ROADMAP.md`.
