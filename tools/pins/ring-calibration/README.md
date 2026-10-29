# tools/pins/ring-calibration: the ring-calibration and q pins

Written before the Swift tests read it; a difference between
the Swift and the pins is a failure of the port, never a band to widen.

```
../References/python/orix/bin/python -I -B tools/pins/ring-calibration/ring-pins.py [out.json]   # about 2 s, deterministic (no clock)
```

- `ring-pins.py` writes `Tests/DiffractoGeometryTests/Fixtures/numpy-ring-calibration-pins.json` (about 143 kB): the CODATA 2018
  constants and the derived hc/e and m0c2/e; the wavelength at 60, 80, 120, 200 and 300 kV; 320 cases (Au {111} {200} {220}
  {311} {422} at a = 4.0782 and Al {111} {200} {220} at a = 4.0495; V, R in {40, 80, 200, 250} px, p in {0.014, 0.05} mm) with
  lambda, 2 theta, L_px, dq and L_mm; 68 refusal plants (each refusal on its boundary, several at once, invalid input); 18
  Bragg-limit plants (lambda/2d = 1/sqrt(2), both sides at 1e-9, 1e-6 and a factor 2, at 60 / 200 / 300 kV); 10 radius plants
  (3 % ellipticity, either axis major); the literals of `RingCalibrationTests` recomputed in decimal (`truth`); the mutation list.
- The oracle is 60-digit `decimal` (arcsin by its power series, tan 2 theta by the identity 2 s sqrt(1 - s^2) / (1 - 2 s^2), so no
  trigonometric library is trusted), with NumPy 2.5.4 routes beside it that must agree to 1e-14. It is written from the
  definitions, ported from macTEM 9482831 `tools/ring-calibration-parity/reference.py` (read by `git show`).
- Exit 0: every self-check holds. 1: a self-check failed. Includes: each mutation (small angle, q = 2/d, Angstrom read as nm,
  pitch in um, non-relativistic lambda, kV as V, V defaulted, radius = a) moves its target quantity by more than 1e-6, so a
  mutation that moves nothing cannot be mistaken for a proof.
- `RingCalibrationPinsTests` (Swift Testing) reads the fixture by `#filePath` and holds every number to 1e-12 relative, every
  refusal by name and in order, `withholdsDeltaQ` as registered. `RingCalibrationTests` are the 30 hand-arithmetic tests of
  macTEM 9482831, converted.
- One literal to know: the typed L_px of Au {111}, R 150 px, 200 kV (14082.307097295722, macTEM's test and
  `RingCalibrationTests`) is one ulp (1.3e-16 relative) from the correctly rounded decimal result (14082.307097295723). It is far inside the
  1e-12 band; the generator holds the typed literals to 1e-15 relative and says so.

## q and its inverse

```
../References/python/orix/bin/python -I -B tools/pins/ring-calibration/q-pins.py [out.json]
```

`q-pins.py`: q and its inverse, 60-digit decimal truth computed two ways (algebraic and series, agreeing to 1e-50) beside
NumPy 2.5.4; writes `Tests/DiffractoGeometryTests/Fixtures/q-pins.json` (byte-identical on rerun). 100 forward, 120 inverse,
320 identity and 51 NaN plants (2026-10-15). It reads `numpy-ring-calibration-pins.json` for the 1/d identity, so
`ring-pins.py` runs first.
