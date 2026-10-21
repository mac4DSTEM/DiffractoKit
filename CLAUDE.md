# DiffractoKit

Crystallography, diffraction geometry, ring and ellipse fits, indexing, orientation and ACOM. Pure Swift, no UI, standalone; the suite's rules are in `../CLAUDE.md`. Start with `README.md`,
`docs/parity.md` and the task.

## Reference
pymatgen 2026.9.24 (spglib 2.8.0, NumPy 2.5.4, SciPy 1.18.1) for cells, symmetry, d-spacings and reflections; py4DSTEM `f050d207` for the scattering factors (bit for bit) and ACOM; orix / diffsims for orientation conventions. Pins and venvs under the gitignored `References/`. Fixtures are generated from it (`tools/*/pins.py`) and committed under `Tests/*/Fixtures/`; a difference inside
the reference's own tolerance is agreement, outside it a bug. Readers, I/O, indices and metadata exact.

## Build and test
```sh
swift build
swift test -c release          # the gate; Debug is many times slower on the numeric loops
tools/consumer-check/run.sh    # a separate package that compiles against the public API
```
Never `swift build --build-tests -c release` (fails on this toolchain). Real-data tests skip when the data is absent.

## Modules
`DiffractoCore` (cells, CIF, symmetry, reflections), `DiffractoScattering` (scattering factors, structure factors), `DiffractoGeometry` (pattern geometry, ring and ellipse fits, calibration, indexing, zone axes), `DiffractoOrientation` (ACOM matcher, orientation results; depends on ImagoKit's `ImagoFFT`).

## Rules of the kit
- Public API only; nothing names an app or an app type; error messages never say "mac4DSTEM" or "macTEM".
- Rows and columns, not x and y. Every accelerated path keeps the CPU float64 path as its reference (D16).
- A new or changed number gets a test against the reference at its upstream tolerance, and a row in `docs/parity.md`.
- Depends on ImagoKit `"0.2.0"..<"1.0.0"` (FFT1D for the ACOM azimuthal transforms). `CrystalPresets` are the preset structures.

## Closeout
`swift test -c release` green; `CHANGELOG.md` entry with a Migration paragraph where a consumer must change something;
version bumped in the entry; one commit on `main` whose last paragraph is the hand-off. Paul tags with
`../macEM-suite/tools/ship.sh` (tag = CHANGELOG top version).
