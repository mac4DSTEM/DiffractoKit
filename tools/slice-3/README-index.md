# tools/slice-3: the indexing pins (RingIndexing, slice 3a)

Registration `docs/registration-slice-3.md` (sections 1, 3, 4, 5). Written before the Swift test reads it; a difference is a
failure of the port, never a tolerance to widen.

```
References/pymatgen-venv/bin/python -I -B tools/slice-3/index-pins.py     # about 2 s, deterministic: writes slice-3-index-pins.json
```

- `index-pins.py` ports macTEM 9482831 `tools/indexing-parity/reference.py` (read by `git show`) onto the KIT's ring lists:
  the families after the |F| cut that slice 2 pinned in `Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json` (so
  `tools/slice-2/pins.py` runs first), at k_max 1.2: Au, Al, Si, Mg, ideal hcp, Weak-cut (macTEM's six) and beta''.
- The decision is NumPy written from the registration text, with a plain-Python scalar re-implementation that must agree on
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
