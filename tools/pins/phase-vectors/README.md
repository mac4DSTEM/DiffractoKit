# tools/pins/phase-vectors: the pin generators

Written before the Swift of slice 4 (DiffractoOrientation) was trusted; the Swift
tests read the file these write. A difference is a failure of the port, never a tolerance to widen (except the changes
named under "Measured floors and open points").

```
../References/python/pymatgen/bin/python -I -B tools/pins/phase-vectors/pins.py     # about 3 s, peak about 330 MB, deterministic
```

- `pins.py` runs on the pymatgen venv (Python 3.14.7, NumPy 2.5.4, SciPy 1.18.1, pymatgen 2026.9.24) and starts
  `pins-py4dstem.py` as a subprocess with `../References/python/py4dstem/bin/python -I -B` (Python 3.12.4, NumPy 2.5.0,
  SciPy 1.18.0, py4DSTEM 0.14.19 = f050d207; read-only, nothing installed into it). The helper checks that the installed
  py4DSTEM files (`crystal.py`, `single_atom_scatter.py`, `utils.py`, the scattering table) are byte equal to the source tree at the
  lock and that the table's sha256 is the lock (`be845d1d...`). One memory cap covers both (the helper is a child). SciPy 1.18.1
  exists only in the pymatgen venv, so the Euler pins are there.
- Output: `Tests/DiffractoOrientationTests/Fixtures/py4dstem-phase-vector-pins.json` (about 1.9 MB, sorted keys, two runs give the same bytes).
- Exit status: 0 all self-checks and predictions hold; 1 a self-check failed (assert); 3 the file was written but a registered or
  predicted number differs from the measured one (printed twice). At the time of writing it exits 3 for the three differences
  listed below; none was edited to fit.

## The three oracles in the file

| Name | What | Where in the JSON |
|---|---|---|
| V2, the reference | NumPy float64 from the written definitions: f_e from the sha-locked py4DSTEM table, F = (1/V) sum f_e occ exp(2 pi i h.x), every hkl with `|F| > 1e-4` (strict) and `|g| <= kMax` enumerated by the bound `|h| <= kMax |a_1|` (cannot miss a lattice point), `s = g.n`, `I = |F|^2 exp(-(s/w)^2)`, kept when `|s| <= slab`, `1e-9 < |q| <= kMax`, `I > 0`; frame `e1 = normalise(n x r)`, `e2 = n x e1`; floor, cap, relative intensity, sort by `|q|` | `phases[].zone_axes[].entries` (pre-cap, pre-floor) and `.library` (after floor and cap) |
| V0, py4DSTEM | `Crystal.generate_diffraction_pattern` at 1e30 V (lambda about 1.2e-26 A), `sigma_excitation_error = w/sqrt(2)`, `tol_excitation_error_mult = slab/sigma`, `tol_intensity = 0`, `k_max = kMax`, `keep_qz`; (0 0 0) dropped; `parse_orientation(zone_axis_lattice=uvw)` | the `py4dstem_*` columns of `entries`, `py4dstem`, `frame_map` |
| V1 | NumPy flat excitation on py4DSTEM's own hkl, g and float32 `|F|^2` | only its residual against V0, `counts.lambda0_vs_flat_numpy_on_py4dstem_table` |
| 200 kV | the same call at `setup_diffraction(200e3)`, curved Ewald sphere | `zone_axes[].kv200` |

`entries` columns (`entry_fields`): h, k, l, gx, gy, gz (the kit's Cartesian g, no 2 pi), s, q1, q2 (the kit frame), |q|, |F|^2,
I_excited, I_relative, then py4DSTEM's x', y', |q|, I_relative (null when py4DSTEM does not list the hkl). `library` columns
(`library_fields`): h, k, l, q1, q2, |q|, I_relative, sorted by |q| (ties: descending hkl; tied |q| are rounding-equal, so a test
matches by hkl).

## What each block is

- `phases`: Al 4.0495 (fcc, Z 13), Au 4.0782, Cu 3.6149 (the apps' presets of mac4DSTEM 7d6fdbd1, the sites as `Crystal.fcc`
  writes them), theta'-Al2Cu, beta''-Mg5Si6 (the apps' `Crystal.betaDoublePrime`, 22 sites in the apps' order and unwrapped
  expressions; `cif_check` shows each site, modulo 1, is a site of the CIF slice 2 pinned, worst difference recorded), Mg hcp
  (3.2094 / 5.2108, sites (0,0,0), (2/3,1/3,1/2)), one triclinic cell. Each carries the cell, the sites with occupancy, `lat_real`,
  `lat_inv`, and per zone axis: `n`, `e1`, `e2`, py4DSTEM's x', y', z', `frame_map` (the 2x2 `M` with `(g.x', g.y') = M (g.e1, g.e2)`,
  `phi_py_rad`, `phi_py_deg`, `det`, `handedness`), `counts`, `flags`, `entries`, `library`, `kv200`.
- Zone axes: cubic Al, Au, Cu [001] [011] [111] [112] [123]; theta' the same plus [100] (edge-on); beta'' the five plus [010]
  (the needle); Mg [001] [100] [110] [210] [111]; triclinic [001] [010] [100] [011] [111]. 37 in all.
- `structure_factor_check` (per phase): py4DSTEM's hkl set and float32 `|F|^2` against V2's, and against slice 2's pins at kMax
  1.2 for Al, Au, Mg, beta''.
- `euler`: 1000 seeded random rotations (SciPy `Rotation.random`, `default_rng(20261015)`), the gimbal and near-gimbal cases
  (Phi = 0, pi, 1e-9, 1e-10, 5e-9, 2e-8, 5e-8, 9e-8, 1e-7, 1e-6, 1e-5, 1e-3, pi - 1e-5, pi - 1e-9, pi - 1e-10, pi - 5e-9, each with three
  (phi1, phi2)) and the 37 py4DSTEM `parse_orientation` matrices (the [001] ones are Phi = 0 exactly). `M` is the py4DSTEM orientation
  matrix, row-major; `R = M^T`; `angles_raw` are SciPy's `as_euler("zxz")`, `angles` wrapped to (-pi, pi]; `crystal_to_lab` is
  `from_euler("zxz", angles).as_matrix()`; `scipy_gimbal_warning` is whether SciPy locked.
- `frame_matrix_cases`: 20 `B R_z(a)` products (10 bases are the columns (e1, e2, n) of frames above, 10 seeded random
  rotations; a in {0, pi/2, -pi, 17 random}).
- `closed_forms`: `chanceMatchFraction` (25 cases), `nearest` (50 seeded + the exact boundary, tie and edge flags), `score` (20 + the
  boundary), `containsFriedelPair` (20 + the boundary), `inPlaneSteps` (13 steps), `lowIndexZoneAxes`, `PhaseVectorResolution`
  (not pinned, see below).

## Measured floors and open points

1. **lambda -> 0 works through py4DSTEM's own API**: `setup_diffraction(1e30)` gives s = g.n exactly in practice. Against NumPy
   flat excitation on py4DSTEM's own table the residual is 0 on six phases and 1.5e-14 relative on beta'' (the `exp(-s^2/(2 sigma^2))`
   versus `exp(-(s/w)^2)` arithmetic), the hkl sets are identical, the projected positions agree to 2.2e-16: below the 1e-12 bar, so
   py4DSTEM's lambda -> 0 pattern is pinned directly (V0) as well as the NumPy value.
2. **py4DSTEM's structure factors are complex64.** `calculate_structure_factors` accumulates `struct_factors` as `complex64`
   (`|F|^2` float32), so its relative intensity differs from the float64 value by up to 1.15e-5 (beta'', the weak reflections),
   3.1e-6 (theta'), 1.8-3.6e-7 on the others, and `|F|^2` by up to 1.3e-5. **G2's "relative intensity within 1e-9" cannot hold against
   py4DSTEM**; adopted instead: hkl sets equal and `|q|` within 1e-12 against py4DSTEM (measured 8.9e-16), relative intensity within
   5e-5 against py4DSTEM (the float32 floor, measured 1.2e-5), within 1e-12 against the float64 reference V2 (the `entries` /
   `library` columns); the kit's q rotated by phi_py equals py4DSTEM's x', y' within 1e-12 (measured 8.9e-16).
3. **|F|^2 against slice 2** (kit float64, Al, Au, Mg, beta'' at kMax 1.2, kept sets equal): max relative difference 1.6e-15 (Al),
   4.3e-16 (Au), 7.3e-16 (Mg), 5.5e-13 (beta'', slice 2 pinned the CIF's wrapped sites). py4DSTEM's float32: 1.7e-7, 1.5e-7, 2.7e-7, 1.1e-5.
   Same convention (Lobato table, Debye-Waller off, occupancy, 1/V, no 2 pi).
4. **The hkl sets**: py4DSTEM's tiling misses no kept hkl on any of the seven phases at kMax 1.6 (beta'' included), and keeps none the
   float64 cut drops; the sets at the zone axes equal V2's in all 37 cases.
5. **phi_py and handedness**: every axis is a rotation (det +1; both frames are right-handed about n). Al, Au, Cu: [001] 90, [011] 0,
   [111] 30, [112] 39.2315, [123] 21.8454 degrees; theta': [001] 90, [100] 180, [011] 0, [111] 35.4662, [112] 41.8949, [123] 23.9296;
   beta'': [001] -90, [011] -19.9543, [111] 54.3291, [112] 64.3455, [123] 45.4439, [010] 0; Mg: [001] 90, [100] 180, [110] 0, [210] 0,
   [111] 26.1782; triclinic: [001] 25.1358, [010] 0, [100] 180, [011] 0.5718, [111] 22.8484. All stored in `frame_map`.
6. **200 kV, curved Ewald sphere** (the same sigma and tolerance): the slab set changes only for theta' [112] (6 hkl enter, 0 leave),
   theta' [123] (21 in, 0 out), beta'' [111] (51 in, 0 out), beta'' [123] (0 in, 54 out); the other 33 zone axes are identical sets
   (relative intensities differ by 1-24%, `max_abs_dI_relative_common`). A parity row, not a gate.
7. **The low-index zone axes: the first prediction said 62, the stated definition gives 49** (98 primitive vectors in the box, divided by
   2 for sign). 62 is the count of the 124 nonzero vectors deduplicated up to sign without the gcd reduction; both lists are in the
   file (`axes`, `alternative_unreduced`). The kit follows the stated definition (49).
8. **SciPy's gimbal rule is not the one first written.** SciPy 1.18.1 locks (sets the third angle to 0) when the middle angle is below
   1e-7 (100% at 9.9e-8, 50% at exactly 1e-7, 0% at 1.01e-7), the written rule was |sin Phi| <= 1e-8. Cases with Phi in {2e-8, 5e-8,
   9e-8, 1e-7} are marked `rules_agree: false`; tests skip them. Phi = 1e-5, 1e-6
   and 1e-3 are unlocked in both, so the threshold mutation (1e-4) is caught.
9. **Definitions not written before the run** (so the pins hold a stated reading, flagged in the file): `chanceMatchFraction` guards
   (R <= 0, non-finite inputs: not pinned), `score` (pinned as the nearest-assignment reading; the any-pair reading is recorded and the
   cases where they differ are flagged; none differ in the 20 seeded cases), `containsFriedelPair` (some pair with |v_i + v_j| <= radius),
   the unit of `inPlaneSteps` rounding (half away from zero). `PhaseVectorResolution` and `scaledToDetector` are **not pinned**: the written
   definitions do not cover them and the oracle may not read the Swift.
10. theta'-Al2Cu has no definition in the apps or the kit. Pinned: body-centred tetragonal 4.04 / 5.80 A, Cu 2a (0,0,0), (1/2,1/2,1/2),
   Al 4d (0,1/2,1/4), (1/2,0,1/4), (0,1/2,3/4), (1/2,0,3/4) (the fluorite-derived arrangement), sites written out. The triclinic cell is
   slice 1's (4, 5, 6, 80, 85, 95) with two invented sites: the kit has no triclinic CIF fixture.
