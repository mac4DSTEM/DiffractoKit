# tools/slice-3b: the ring and ellipse fit (DiffractoGeometry, slice 3b)

Registration: `docs/registration-slice-3.md` sections 2, 3 and 5. Seeds read only by `git show`: mac4DSTEM `7d6fdbd1`
(`EllipseCalibration.swift`, `DiffractionPattern.swift`, `Calibration.swift`, the noisy, stop-rule and row-29 tests, the
`tools/ellipse-calibration-test` sweep) and macTEM `9482831` (`RingCalibration+Fit.swift`, its 8 tests).

- `apps-parity/run.sh` (+ `main.swift`): compiles the apps' own `EllipseCalibration.swift` and `DiffractionPattern.swift` at 7d6fdbd1
  with `swiftc -O`, renders 35 planted rings (3 + 4 noisy, the 4 stop-rule rings, 20 clean sweep rings, 4 sparse rings) with the
  SplitMix64 + Box-Muller renderer, runs `fit1D` (three ways), `fitAmorphousRing` and `fitBestAvailable` (two ways) on each and
  writes every field as bit patterns to `Tests/DiffractoGeometryTests/Fixtures/slice-3b-apps.json`; the sweep rings' raw Float32
  files go to `.build/apps-parity-3b/sweep`. The test `EllipseAppsParityTests` requires the kit equal in every bit.
- `sweep.py` (`../References/python/orix/bin/python -I -B tools/slice-3b/sweep.py`): SciPy 1.18.1 `least_squares` optimum of fit1D's
  objective on the 20 sweep rings -> `Fixtures/slice-3b-sweep.json`; held in `EllipseCalibrationTests`.
- `mutants.py <copy-dir> <build-dir>`: the break-once run, in a copy of the repository.

Order: `apps-parity/run.sh` first (it writes the dump and the pixel files), then `sweep.py`.

## Results (2026-10-15, DiffractoKit lane 3b)

- Bit identity: the kit's `fit1D` / `fitAmorphousRing` / `fitBestAvailable` equal the apps' (mac4DSTEM 7d6fdbd1) in every output bit
  on 35 rings x 6 calls, compiled -O and compiled -Onone (run.sh kit-check, and `EllipseAppsParityTests` in a Debug `swift test`
  against the -Onone dump; a release `swift test` reads the -O dump). The apps' own -O and -Onone dumps differ in the last bits of
  some profile fits (25 ulps in a centre, 98 in a theta: reassociation/contraction at -O), so each build is held to the apps built
  the same way; this is why there are two fixtures.
- The parity set does not see three of the nine mutants (floor 0.30, the profile residual cap, the outer annulus boundary: no ring
  of the set ends within 0.30 px of the floor, or at the cap, or has a pixel exactly on the boundary); `EllipseCalibrationTests`
  holds those three (`parity-mutants.py` shows which).
- macTEM's 8 `RingFitCalibrationTests` pass unchanged on the 7d6fdbd1 fit. Against macTEM's pre-ADR-069 copy the planted ring's fit
  moved by 1e-9 px / 1e-9 relative (clean 2 %: a 60.594059408098644 against 60.59405940811859) and the 12 % ring likewise.
- Clean sweep against SciPy 1.18.1: centre 7.1e-8 px, axes 1.2e-9 relative, theta 5.0e-8 rad (bands 1e-5, 1e-7, 1e-6).
