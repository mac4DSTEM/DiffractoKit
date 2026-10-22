# Status

DiffractoKit 0.10.0 (2026-10-22), on ImagoKit `"0.2.0"..<"1.0.0"`. Live: cells (seven systems, Niggli, conventional cell), CIF import
with symmetry detection and a completeness check by IT number or by symbol, space-group-free symmetry from the atoms, reflections and rings, kinematic
structure factors (py4DSTEM table bit for bit), pattern geometry (ring and ellipse fits, calibration, zone axes, indexing), the
ACOM matcher (CPU reference, Metal option, optional in-plane refinement), vector-matched phase maps, orientation symmetry
and IPF keys, misorientation between two groups, orientation-relationship variants (beta'' in Al: 12), orientation-map planes,
the ring-sequence assignment of the known-crystal Q, delta q from several rings, and single-pattern beam centre, radial profile and spot refinement. What is checked against what: `docs/parity.md`; timings against Python: `docs/benchmarks.md`;
what is open: `docs/open-items.md`.

Tests are split in two (ROADMAP task 31): `DiffractoKitParityTests` holds the pinned comparisons against the Python references (slice pins for cells and structures, scattering, ring calibration and indexing, ellipse apps parity, slice 4 pins and apps parity, py4DSTEM family and orix pins; 91 tests), the other targets hold the fast unit tests (519 tests), and `DiffractoKitTestSupport` holds the helpers both use (`FixtureDirs`, the pin decoders, structure builders; the fixtures stay in the Fixtures folders of the older test targets). Measured on the Mac (2026-10-21, warm build): `swift test -c release --parallel --skip DiffractoKitParityTests` 15 s, the full `swift test -c release --parallel` 21 s, before the split 19 s; 610 tests before and 519 + 91 = 610 after.
