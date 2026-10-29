# tools/pins/ellipse: the ring and ellipse fit (DiffractoGeometry)

Seeds, from the apps' history before they were public: mac4DSTEM `7d6fdbd1`
(`EllipseCalibration.swift`, `DiffractionPattern.swift`, `Calibration.swift`, the noisy, stop-rule and row-29 tests, the
`tools/ellipse-calibration-test` sweep) and macTEM `9482831` (`RingCalibration+Fit.swift`, its 8 tests).

- `sweep.py` (`../References/python/pymatgen/bin/python -I -B tools/pins/ellipse/sweep.py`): renders the 20 clean sweep rings (a Python
  port of the planted ring of `Tests/DiffractoKitTestSupport/EllipseFixtures.swift`, hash-equal to the Swift renderer) and writes
  SciPy 1.18.1's `least_squares` optimum of fit1D's objective on each, with its plant and start, to `Fixtures/scipy-ellipse-sweep-pins.json`;
  held in `EllipseCalibrationTests`.

## Results (2026-10-15)

- Bit identity (2026-10-15, a one-time port check; the dump and its test were retired in #18): the kit's `fit1D` /
  `fitAmorphousRing` / `fitBestAvailable` equalled the apps' (mac4DSTEM 7d6fdbd1) in every output bit on 35 rings x 6 calls,
  compiled -O and -Onone.
- That parity set did not see three of the nine mutations (floor 0.30, the profile residual cap, the outer annulus boundary: no ring
  of the set ends within 0.30 px of the floor, or at the cap, or has a pixel exactly on the boundary); `EllipseCalibrationTests`
  holds those three.
- macTEM's 8 `RingFitCalibrationTests` pass unchanged on the 7d6fdbd1 fit. Against macTEM's pre-ADR-069 copy the planted ring's fit
  moved by 1e-9 px / 1e-9 relative (clean 2 %: a 60.594059408098644 against 60.59405940811859) and the 12 % ring likewise.
- Clean sweep against SciPy 1.18.1: centre 7.1e-8 px, axes 1.2e-9 relative, theta 5.0e-8 rad (bands 1e-5, 1e-7, 1e-6).
