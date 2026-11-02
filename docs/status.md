# Status

DiffractoKit 0.10.0 (2026-10-22), on ImagoKit `"0.2.0"..<"1.0.0"`.

Live:

- Cells in the seven systems, Niggli reduction, the conventional cell.
- CIF import with symmetry detection and a completeness check by IT number or by symbol; space-group-free symmetry from the
  atoms; reflections and rings.
- Kinematic structure factors (py4DSTEM's scattering-factor table, bit for bit).
- Pattern geometry: ring and ellipse fits, calibration, zone axes, indexing; the ring-sequence assignment of a known crystal,
  Δq from several rings, and the single-pattern beam centre, radial profile and spot refinement.
- The ACOM matcher (CPU reference, Metal option, optional in-plane refinement) and vector-matched phase maps.
- Orientation symmetry and IPF keys, misorientation between two groups, orientation-relationship variants (beta'' in Al: 12),
  orientation-map planes.

Tests: `DiffractoKitParityTests` holds the pinned comparisons against the Python references; the other test targets hold the
fast unit tests; `DiffractoKitTestSupport` holds the helpers both use.

What is checked against what: `docs/parity.md`. Timings against Python: `docs/benchmarks.md`. What is open: the GitHub
issues; known limits: `docs/parity.md`. Where it is going: `ROADMAP.md`.
