# tools/slice-3: the ring-calibration and indexing pins

Written before the Swift tests read it; a difference between
the Swift and the pins is a failure of the port, never a band to widen.

```
../References/python/orix/bin/python -I -B tools/slice-3/ring-pins.py [out.json]   # about 2 s, deterministic (no clock)
```

- `ring-pins.py` writes `Tests/DiffractoGeometryTests/Fixtures/slice-3-ring-pins.json` (about 143 kB): the CODATA 2018
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

## The indexing pins

Written before the Swift test reads it; a difference is a
failure of the port, never a tolerance to widen.

```
../References/python/orix/bin/python -I -B tools/slice-3/index-pins.py     # about 2 s, deterministic: writes slice-3-index-pins.json
```

- `index-pins.py` ports macTEM 9482831 `tools/indexing-parity/reference.py` (read by `git show`) onto the KIT's ring lists:
  the families after the |F| cut that slice 2 pinned in `Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json` (so
  `tools/slice-2/pins.py` runs first), at k_max 1.2: Au, Al, Si, Mg, ideal hcp, Weak-cut (macTEM's six) and beta''.
- The decision is NumPy written from the definition in `RingIndexing.swift`'s header, with a plain-Python scalar re-implementation that must agree on
  every plant (z = |q - g|/sigma, delta = g/q - 1, candidates by (z, g), within = z <= k, ambiguous on two within or a margin
  below ln 19, reach = sqrt(k^2 + 2 ln 19) sigma, the verdict order invalidInput, ambiguous, beyondKMax, none, accepted).
- **T1**: every family x u in {1e-4, 1e-3, 3e-3, 1e-2, 3e-2} x s in {+-2.9, +-2, +-1, +-0.5, 0}; q = g(1 + s u), sigma = u q,
  k = 3 (22,500 plants). Closed forms P1 are checked in the generator (z = |s|/(1 + s u), delta = -s u/(1 + s u)).
- **T4**: 35 refusal plants on the kit's lists (invalid numbers, beyondKMax, the one-family list Au at k_max 0.45 and the
  empty one at 0.3, forbidden positions, ambiguous pairs, the margin plant, beta'' pair and its accepted {2 0 0}).
- Self-checks (exit 1 on failure, the file still written): the T1 counts over macTEM's six lists equal its registered
  prediction per u (711/0/0/36/0 ... 131/19/572/11/14); beta'' accepted counts per u equal 2584, 508, 184, 58, 25; beta'' at
  s = 0 has 23, 8 and 4 accepted families at u = 3e-3, 1e-2 and 3e-2 (macTEM ADR 024) and the eight names at 1e-2; every T4
  plant carries the verdict it was planted for; the closest plant to z = k, margin = ln 19 and q + R = k_max is 6.7e-6, 6.2e-4
  and 2.9e-4 away (macTEM's redraw distances 1e-6, 1e-6, 1e-9).
- Read by `Tests/DiffractoGeometryTests/RingIndexingPinsTests.swift`: verdict and payload family names equal, delta 1e-14,
  z 1e-10 max(1, z), margin 1e-9 max(1, margin), reach 1e-9 max(1, reach).
- Non-finite numbers are written as the strings "NaN", "Infinity", "-Infinity".
