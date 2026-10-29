# Orientation symmetry and colour conventions

`CrystalOrientationSymmetry(laueClass:)` supports every `LaueClass`. Its operators are proper rotations in a standard
Cartesian frame: monoclinic b is y; tetragonal c is z; trigonal/hexagonal a is x, c is z and b is at 120° to a. The
trigonal -3m has two settings in hexagonal axes, chosen with `CrystalOrientationSymmetry(laueClass: .threeBarM, setting:)`
(`ThreeBarMSetting`): `.threeBarOneM` (-31m: P-31m, P-31c; twofolds perpendicular to a, at -30, 30 and 90 degrees from x; pymatgen's
`PointGroup("-3m")`; the default) and `.threeBarMOne` (-3m1: P-3m1, P-3c1, R-3m and R-3c in hexagonal axes; twofolds on a, b and -(a+b),
at 0, 60 and 120 degrees; orix's `D3d`). A cell in another setting must first be transformed into one of these frames
(a rhombohedral-axes R-3m cell is not covered). This constructor does not infer a setting from a cell or from its atoms;
`PhaseModel.orientationSymmetry(of:)` does.

| Laue class | Proper subgroup order | Generators (fractional for hexagonal axes) |
|---|---:|---|
| -1 | 1 | identity |
| 2/m | 2 | C2(y) |
| mmm | 4 | C2(x), C2(y) |
| 4/m | 4 | C4(z) |
| 4/mmm | 8 | C4(z), C2(x) |
| -3 | 3 | C3(z) |
| -3m | 6 | C3(z), pymatgen's basal twofold |
| 6/m | 6 | C6(z) |
| 6/mmm | 12 | C6(z), C2(a) |
| m-3 | 12 | C2(x), C3([111]) |
| m-3m | 24 | C4(z), C3([111]) |

`reduce(_:)` chooses the minimum-angle `S * matrix`, where the input's columns are lab axes in crystal coordinates
(the py4DSTEM matrix convention). The returned matrix is one proper rotation. Canonical quaternion `(−w,x,y,z)`
comparison breaks ties, ignoring component differences up to 1e-13. Near a half-turn (`|w| < 1e-14`), the first imaginary
component whose magnitude exceeds 1e-14 fixes the sign. Otherwise w is nonnegative. This is the minimum-angle
orix/MTEX convention; orix 0.15.0 is executed as an oracle (MTEX is MATLAB and is not). The fixture independently minimizes over
pymatgen's operation sets with SciPy quaternions, and (0.6.0) orix 0.15.0's `Orientation.reduce()` and `angle_with` agree with it for
all eleven classes (see "Disorientation" below). The diagnostic angle is Float radians.

```swift
import DiffractoCore
import DiffractoOrientation

let symmetry = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
let reduced = try symmetry.reduce(euler.py4DSTEMOrientationMatrix)
let beams = try symmetry.sampleFundamentalZone(count: 64)
let nativeRGB = try symmetry.ipfColor(direction: reduced.matrix.columns.2, key: .native)
```

`reduceDirection` folds poles modulo the proper operators and inversion. The original native cubic, hexagonal and
monoclinic sectors are preserved. Other classes choose the largest `(z,x,y)` orbit member with component comparisons
at 1e-13. Sampling uses folded Fibonacci candidates and spherical farthest-point selection; the three original native
samplers preserve their seeds, candidate counts and arithmetic. Other classes seed folded z, x, y (deduplicated), then
use at least 192 candidates or 16 per requested point. This is a deterministic coverage policy, not py4DSTEM's angular
sampling grid or a claim of optimal covering radius. Counts outside 0...4096 throw before allocation.

Both are cross-checked against orix 0.15.0 (`tools/pins/orix/direction-pins.py`, `OrixDirectionTests`): the kit's zone is not orix's (orix's
`fundamental_sector` is a different representative region), so the comparison is of orbits under orix's own Laue group in the kit frame. On
3 765 poles (300 seeded random and 20 designed per class, twelve class/setting entries) the reduced pole lies in the orbit of the input and of orix's
folded pole (worst gap 1.1e-15), and every member of an orbit reduces to the same pole (9e-16), so the kit's zone holds exactly one pole per
orbit. Sampled poles lie in the zone, no two are equivalent under orix's group, and the covering radius over 1 500 Fibonacci probe poles
(radians, counts 1 / 3 / 8 / 16 / 64) is: -1 1.570 / 0.950 / 0.858 / 0.499 / 0.266; 2/m 1.571 / 0.950 / 0.587 / 0.368 / 0.187; mmm 1.570 / 0.950 / 0.473 / 0.335 / 0.133;
4/m 1.570 / 0.669 / 0.465 / 0.319 / 0.134; 4/mmm 1.570 / 0.656 / 0.320 / 0.241 / 0.096; -3 1.570 / 0.800 / 0.418 / 0.312 / 0.152; -31m 1.570 / 0.800 / 0.380 / 0.256 / 0.104;
-3m1 1.570 / 0.800 / 0.339 / 0.229 / 0.110; 6/m 1.570 / 0.800 / 0.321 / 0.232 / 0.106; 6/mmm 1.570 / 0.800 / 0.266 / 0.191 / 0.083; m-3 0.950 / 0.476 / 0.307 / 0.169 / 0.080;
m-3m 0.950 / 0.480 / 0.183 / 0.128 / 0.059 (a deterministic coverage policy, not an optimum). No change was needed.

## Three keys with different input contracts

All keys return RGB as `SIMD3<Float>`, without UI or a map container. The convention must be selected explicitly.

- `.native` accepts the original Cartesian pole and folds it. For m-3m (TSL-style equalized spherical key), 6/mmm
  (square-root triangular key) and b-unique 2/m (HSV disk) it preserves the app's current **own convention**, pinned to the
  app's Float bits. For the other eight classes (0.6.0) it is the orix key below, so every class has a native key and
  `nativeKeyUnavailable` is no longer thrown (the case remains so exhaustive `switch` statements compile).
- `.orix` (0.6.0) accepts the original Cartesian pole and colours it with orix 0.15.0's `DirectionColorKeyTSL` (the EDAX/TSL
  key as implemented in MTEX) for **all eleven** classes: the pole is folded into orix's `fundamental_sector` with
  `Vector3d.in_fundamental_sector` (MTEX `project2FundamentalRegion`), then coloured by the polar coordinates about the
  sector centre (the 1000-point corrected azimuth for three-vertex sectors, hue from azimuth, lightness `0.5 + polar/2`).
  The sectors, the centres, orix's element order for ties and the frame below are a generated table
  (`OrixIPFSectors.swift`, `tools/pins/orix/ipf-pins.py`). Measured against orix on 2 777 poles (240 seeded signed random poles
  per class plus the axes, sector vertices and centres, and points beside the vertices): worst RGB gap 2.98e-8, the Float
  rounding (the test tolerance is 1e-7). The measured distance between the app's native keys and orix's is recorded by a
  test: m-3m 3.5e-6 (the app's key is the same construction), 6/mmm 0.30 and 2/m 1.0 (different conventions).
- `.py4DSTEM` accepts a pole already in py4DSTEM's `Orientation.family` for the class's **automatic plan sector**.
  The exception is 2/m: pass the **original matrix pole**, because `plot_orientation_maps` uses NNLS there. This key
  matches the full-mask (`mask = 1`) beam RGB from the pinned routine, converted to Float. No pole normalization is
  inserted: the Python clipping and maximum normalization are preserved. Correlation masks, rotated in-plane colours,
  custom plan sectors and plotting are not implemented. Negative sector weights clip to zero, as in Python.

### orix frames

orix has one setting per group and the kit's standard settings differ for two classes. The generator maps orix's group onto
the kit's operators (asserted as sets for all eleven) with a rotation `frame`; the key is orix's, evaluated in the kit's
setting:

| Laue class | orix setting | `frame` (orix to kit) |
|---|---|---|
| 2/m | twofold on z (c unique) | cyclic: kit (x, y, z) = orix (Y, Z, X), so the kit's b axis is orix's z |
| -3m (-31m) | twofolds on x, 60 degrees apart | Rz(-30 degrees): the kit's pymatgen setting has its twofolds at -30, 30 and 90 degrees |
| -3m (-3m1) | the same | identity: the twofolds already lie on a, b and -(a+b) |
| all others | the same | identity |

Two frames are equally valid for -31m (Rz(-30) and Rz(+30) differ by Rz(60), which is not in the group) and for 2/m (any
rotation about b); they give the same colour wheel mirrored. The choices above are fixed and pinned; a consumer that uses
another setting must transform the pole first.

## Disorientation

`CrystalOrientationSymmetry.disorientation(from:to:)` (0.6.0) is the smallest rotation angle between two orientations of one
crystal over their symmetry equivalents: with the py4DSTEM matrix convention of `reduce(_:)` (which is also orix's
`Orientation` matrix: lab to crystal, symmetry on the left) the misorientation is `to * from^T` and the result is the
minimum-angle `S * to * from^T`. Because `angle(S1 M S2)` equals `angle(M S2 S1)`, the two-sided reduction of orix's
`Misorientation((sym, sym)).reduce()` and the one-sided one are the same angle. The cross-check against orix 0.15.0
(`tools/pins/orix/disorientation-pins.py`; 150 seeded random orientations, 10 designed ones (the identity and small, medium and large rotations) and 150 pairs per class, all eleven
classes, `OrixDisorientationTests`): `reduce`'s Float angle within 1.19e-7 rad, the reduced matrix within 6.7e-16 wherever the
minimum is unique (next-best angle more than 1e-6 rad away), the Double pair angle within 1.4e-14 rad. The kit's reduction
needed no change. Notes from the comparison:

- orix 0.15's `Quaternion.from_matrix` zeroes a small scalar part near a half-turn (an error of 1.3e-5 rad at
  |w| = 7e-6); the pins take the quaternion from SciPy instead and require orix's own `to_matrix` to return the matrix.
- The angle of a single orientation to the reference is not invariant under a relabelling of the crystal axes, so for 2/m and
  -3m the pins relabel the lab axes together with the crystal's (`M_orix = frame^T M frame`); pairs need no such step.

## py4DSTEM direction families

`py4DSTEMDirectionFamily(_:order:)` (0.6.0) is py4DSTEM f050d207's `Crystal.symmetry_reduce_directions`, returning `Orientation.family`:
each column of the matrix is folded independently by the first symmetry operator that puts it inside the class's automatic
orientation-plan sector (all three barycentric coordinates `>= 0`), so the result is **not** a rotation and is not the
minimum-angle matrix. Quirks are kept: a column no operator maps into the sector keeps the first operator (-1, 2/m and -3m, whose
automatic sector is not a fundamental domain, return it for part of the sphere), the first qualifying operator wins, and -3
tests its operators twice, the second copy against the sector with its first vertex negated (py4DSTEM's 180 degree fibre
case). The operator list is selected with `order:`. `.pymatgen` (the default) is a real py4DSTEM crystal's: `Crystal.orientation_plan` takes the full Laue group from
`SpacegroupAnalyzer(structure).get_point_group_operations()` (pymatgen's order, identity first, duplicates removed; transposed to
`lat_inv.T @ R.T @ lat_real`), twice over for -3. The kit ships the order of the first space group of each class (P-1, P2/m, Pmmm, P4/m,
P4/mmm, P-3, P-31m for -31m, P-3m1 for -3m1, P6/m, P6/mmm, Pm-3, Pm-3m) as a generated table (`PymatgenOperatorTable.swift`,
`tools/pins/orientation-symmetry/family-real-pins.py`); other space groups of a class that list the operators differently (Fddd; the origin-choice-2 groups of 4/m, 4/mmm,
m-3, m-3m) gave identical families on all pinned matrices. `.sortedFractional` is the pin harness's: the proper operators sorted by their fractional
matrices, then their negatives. They differ in the fallback operator, which is the first one: for 2/m and -3m (-31m) a column that no
operator maps into the automatic sector (which is not a fundamental domain there) keeps the first, a twofold in the sorted list and the identity in a real crystal's.
Pinned on real crystals (57 `Crystal`s of centrosymmetric space groups, 8 607 matrices, `Py4DSTEMRealFamilyTests`): family entries within 4.4e-16, all twelve class/setting
entries; against the same crystals the sorted order disagrees on 139 of 151 matrices for 2/m and 136 of 151 for -31m. Under the harness order, pinned
(`tools/pins/orientation-symmetry/family-pins.py` runs the real function; `Py4DSTEMFamilyTests`): 2 387
matrices, all eleven classes, 7 161 columns, family entries identical (worst gap 0.0); 171 columns lying within 1e-12 of a
sector wall (the identity's columns on the 6/m wall, for example) are not compared, because whether the exact `>= 0` test
passes there is BLAS rounding noise in NumPy. For example, native cubic folding uses [001]-[101]-[111], while the automatic
Python sector uses [001]-[011]-[111].

## Orientation class of a structure

`PhaseModel.orientationClass(of:)` (0.6.0, `DiffractoOrientation`) names the Laue class the reduction can use for a structure:
the class detected from the atoms, when the structure's own proper Laue operators in its cell basis are exactly the class's
standard operators in one of its settings. It agrees with spglib's Laue class for all 15 pinned space groups (-1 to m-3m). For -3m the
structure's setting is named too: `PhaseModel.orientationSymmetry(of:)` returns the `CrystalOrientationSymmetry` with the setting
(-31m for P-31m and P-31c, -3m1 for P-3m1, P-3c1, R-3m and R-3c in hexagonal axes, pinned to pymatgen / spglib's operations of space groups 162 to 167 in
`TrigonalSettingTests`); a rhombohedral-axes cell or a primitive cell of a centred lattice still returns nil.

The py4DSTEM key's automatic-sector inverse matrices are a small generated table from the executed pinned plan parser.
The monoclinic default sector is [001], [010], [100]; NNLS of its positive/negative Cartesian basis is exactly the positive
and negative parts of the pole, so no iterative solver or new dependency is needed. Fused multiply-add retains NumPy's
boundary residual in the other keys. Exactness means Float bit patterns after converting the Python Double result,
not bit equality between different floating-point types. The measured Double-to-Float gap is reported in parity.md.

## Reproducing the fixtures

From the worktree, with the existing pinned references (no installation or network access):

```sh
../../References/python/pymatgen/bin/python -I -B tools/pins/orientation-symmetry/pins.py
swift test -c release --jobs 4
```

The orix pins (`tools/pins/orix/ipf-pins.py`, `tools/pins/orix/disorientation-pins.py`) use the orix venv (orix 0.15.0, NumPy 2.5.4, SciPy 1.18.1, Python 3.12, freeze file in
`References/python/orix/freeze.txt`; see `CONTRIBUTING.md`); the direction-family pins use the py4DSTEM venv. `MACEM_SUITE_ROOT` can locate references from a differently placed checkout. `pins.py` invokes the pinned py4DSTEM
interpreter and compiles only a `git show` extraction of the app's symmetry classes at c809679c into this worktree's
`.build/symmetry-generator`. It does not build or modify the app. The JSON contains seed 20261016, package versions,
source digests (including pymatgen symm_data.json), every proper group, random rotations, half-turn/boundary cases, original native Float colour bits and
sample positions, executed py4DSTEM RGB and Float bits, and the new sectors/samples computed from pymatgen orbits.
Python's family reducer is called only by the generator to obtain representative colour inputs.

The generator enforces NumPy 2.5.4 / SciPy 1.18.1 in the pymatgen environment and NumPy 2.5.4 / SciPy 1.18.1 in the py4DSTEM environment. It also regenerates the colour-basis Swift table.
