# Orientation symmetry and colour conventions

`CrystalOrientationSymmetry(laueClass:)` supports every `LaueClass`. Its operators are proper rotations in a standard
Cartesian frame: monoclinic b is y; tetragonal c is z; trigonal/hexagonal a is x, c is z and b is at 120° to a. The
trigonal -3m setting is pymatgen's standard hexagonal setting. A cell in a different setting must first be transformed
into this frame. This constructor does not infer a setting from a cell or from its atoms.

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
orix/MTEX convention; orix/MTEX themselves are not executed as oracles here. The fixture independently minimizes over
pymatgen's operation sets with SciPy quaternions. The diagnostic angle is Float radians.

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

## Two keys with different input contracts

Both keys return RGB as `SIMD3<Float>`, without UI or a map container. The convention must be selected explicitly.

- `.native` accepts the original Cartesian pole and folds it. It preserves the app's current **own convention** for
  m-3m (TSL-style equalized spherical key), 6/mmm (square-root triangular key), and b-unique 2/m (HSV disk). It does
  not claim numerical parity with orix or py4DSTEM. The other eight classes throw `nativeKeyUnavailable(laueClass)`.
- `.py4DSTEM` accepts a pole already in py4DSTEM's `Orientation.family` for the class's **automatic plan sector**.
  The exception is 2/m: pass the **original matrix pole**, because `plot_orientation_maps` uses NNLS there. This key
  matches the full-mask (`mask = 1`) beam RGB from the pinned routine, converted to Float. No pole normalization is
  inserted: the Python clipping and maximum normalization are preserved. Correlation masks, rotated in-plane colours,
  custom plan sectors and plotting are not implemented. Negative sector weights clip to zero, as in Python.

py4DSTEM's `symmetry_reduce_directions` is **not ported (returns families, not a rotation)**. It transforms the columns
independently; its `family` array need not be orthonormal. Neither the minimum-angle matrix nor `reduceDirection` should
be substituted for that array when reproducing a Python map. For example, native cubic folding uses [001]-[101]-[111],
while the automatic Python sector uses [001]-[011]-[111].

The py4DSTEM key's automatic-sector inverse matrices are a small generated table from the executed pinned plan parser.
The monoclinic default sector is [001], [010], [100]; NNLS of its positive/negative Cartesian basis is exactly the positive
and negative parts of the pole, so no iterative solver or new dependency is needed. Fused multiply-add retains NumPy's
boundary residual in the other keys. Exactness means Float bit patterns after converting the Python Double result,
not bit equality between different floating-point types. The measured Double-to-Float gap is reported in parity.md.

## Reproducing the fixtures

From the worktree, with the existing pinned references (no installation or network access):

```sh
../../tools/memcap.sh -- ../../References/DiffractoKit/pymatgen-venv/bin/python -I -B tools/symmetry/pins.py
swift test -c release --jobs 4
```

`MACEM_SUITE_ROOT` can locate references from a differently placed checkout. `pins.py` invokes the pinned py4DSTEM
interpreter and compiles only a `git show` extraction of the app's symmetry classes at c809679c into this worktree's
`.build/symmetry-generator`. It does not build or modify the app. The JSON contains seed 20261016, package versions,
source digests (including pymatgen symm_data.json), every proper group, random rotations, half-turn/boundary cases, original native Float colour bits and
sample positions, executed py4DSTEM RGB and Float bits, and the new sectors/samples computed from pymatgen orbits.
Python's family reducer is called only by the generator to obtain representative colour inputs.

The generator enforces NumPy 2.5.4 / SciPy 1.18.1 in the pymatgen environment and NumPy 2.5.0 / SciPy 1.18.0 in the py4DSTEM environment. It also regenerates the colour-basis Swift table.
