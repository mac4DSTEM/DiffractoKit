# DiffractoKit 0.3.0: what mac4DSTEM and macTEM change to re-pin

Written 2026-10-15 by DiffractoKit's session, for mac4DSTEM's and macTEM's sessions. 0.3.0 is slice 4: vector-matched phase
mapping and the value types of an orientation result move into `DiffractoOrientation` (registered in
`docs/registration-slice-4.md`, with its dated amendments). It contains 0.2.0 (read `docs/consumers-0.2.0.md` first: the
policies, the floor macOS 15 / iOS 18, the refuter's fixes). **TAG COMMIT: the commit after slice 4's Gate B refuter closes;
filled at the tag.** The tag is the owner's; until it exists nothing here can be pinned.

## In one table

| Who | Deletes | Replaced by |
|---|---|---|
| mac4DSTEM | `mac4DSTEM/Core/Crystal/PhaseVectorMatching.swift` | `DiffractoOrientation`: `PhaseVectorSettings`, `PhaseVectorResolution`, `PhaseVerdict`, `PhaseVectorResult`, `PhaseMap`, `PhaseVectorMatcher`, `PhaseDiskClaim`, `PeakCoordinates` |
| mac4DSTEM | `mac4DSTEM/Core/Crystal/PhaseReferenceLibrary.swift` | `ReferenceVector`, `PhaseOrientationReference`, `LatticeVector`, `OrientationRelationship`, `PhaseDefinition`, `PhaseReferenceSettings`, `PhaseReferenceLibrary` |
| mac4DSTEM | from `mac4DSTEM/Core/Analysis/OrientationResult.swift`: `EulerAngles`, `ACOMOrientation`, the `OrientationResult` struct (one position) | `EulerAngles`, `ZoneAxisFrame` (the renamed `ACOMOrientation`: `detectorBasis(zoneAxis:)`, `matrix(detectorBasis:inPlaneAngle:)`), `OrientationResult` |
| mac4DSTEM | tests: `PhaseVectorMatchingTests`, `PhaseVectorGateBTests`, `PhaseVectorResolutionTests` (classes in `mac4DSTEMTests/PhaseVectorMatchingTests.swift`), `PhaseReferenceLibrarySlabTests.swift` | the kit's `Tests/DiffractoOrientationTests` (57 + 2 of them ported with every expectation unchanged) |
| mac4DSTEM | keeps | `PhaseMapPresentationTests` and the zone-axis / relationship text parsers (app types), `OrientationPlan`, `OrientationMatcher` (template matching: slice 5), the map and image presentation of `OrientationResult.swift`, the cubic symmetry and IPF key |
| macTEM | nothing | it may use `PhaseReferenceLibrary` and `PhaseVectorMatcher.classify` on one pattern's vectors, and `EulerAngles` / `ZoneAxisFrame` for orientation read-outs, instead of writing its own |

## API changes for the switch

- **`Crystal` → `CrystalStructure`.** `PhaseDefinition.crystal` is a `CrystalStructure`. The library reads
  `PhaseReferenceLibrary.reflections(of:kMaxInvAngstrom:)`: the kit's kinematic list at the default |F| cut, re-ordered to the
  app's order (|g|, then h, k, l ascending), bit identical to `Crystal.reflections(kMax:)`. `unsupportedElements(of:)` replaces
  `Crystal.unsupportedElements`. A structure the kit refuses (a close contact, an invalid cell or occupancy, a refused kMax)
  throws `PhaseReferenceLibrary.Failure.structureRefused(phase:reason:)`.
- **`BraggVectors` → `scan: [[SIMD2<Double>]], scanWidth:, scanHeight:`** (row-major, one array per position).
  `PeakCoordinates.detectorPixels(originX:originY:invAngstromPerPixel:)` reproduces the app's `Float` subtraction exactly (pass
  each `BraggPeak`'s `x`, `y` as `SIMD2<Double>(Double(x), Double(y))`); `.calibrated` takes vectors already in Å⁻¹ about the
  origin. The overloads with the app's `originX:originY:invAngstromPerPixel:` labels remain on `map`, `fitMatrixOrientation`,
  `fitZoneAxis`, `claims` and `experimentalVector(s)`.
- **`AnalysisCancellationToken` → closures**: `progress: (@Sendable (Double) -> Void)?`, `isCancelled: (@Sendable () -> Bool)?`
  (`build`, `map`, the fits). `map` gained `maximumWorkers: Int? = nil` (nil = every core, as before; the result does not
  depend on it).
- **`fitZoneAxis` throws** `PhaseReferenceLibrary.Failure` (it builds reflections through the kit).
- Both settings structs have `public static let default`; the 24 values are the app's (`PhaseVectorSettings` 17,
  `PhaseReferenceSettings` 7), each documented with the measurement that set it.

## What numbers move

None on any input the app's code handles deterministically: the kit equals mac4DSTEM 7d6fdbd1's compiled code bit for bit
on about 1.68 million values (gate G1: libraries, 72 maps, fits, claims, 270 Euler cases). Four documented differences:
1. **A rotation about the beam with R22 rounding to 1 − 2⁻⁵³** (and R20 = R21 = 0) was lost by the app (both Euler angles 0,
   the matrix 0.479 off); the kit folds it into φ1 as at R22 = 1. Rare in practice; a fix.
2. **An exact score tie between two phases** was broken by the process's hash seed in the app (the same pattern could be
   labelled differently in two runs); the kit gives it to the lower phase index.
3. **A malformed scan** (position count ≠ width × height) gave a matrix fit in the app's `fitMatrixOrientation` and
   `fitZoneAxis`; the kit returns nil / `[]`, as `map` always did. An out-of-range `matrixEntryIndex` makes `map` return nil
   where the app trapped.
4. **The low-index axis set has 49 directions** (the app's own comment said 62; the code always made 49).

## Checks behind it

G1 apps parity (bit for bit, -O); G2 py4DSTEM f050d207 (hkl sets, |q| 1.1e-15, the frame rotation per axis, intensity within
py4DSTEM's complex64 floor 1.15e-5); G3 SciPy 1.18.1 Euler on 1000 rotations and the gimbal cases; G4 closed forms; G5 the
app's planted-truth tests; every new test broken once; the registered mutants (`tools/slice-4/mutants.py`); one Gate B
refuter on the frozen commit (results filled at the tag). Not validated on a real pattern beyond what the app measured; the
method is Thronsen et al., Ultramicroscopy 255 (2024) 113861 (CC BY 4.0), the paper only.

## Not in 0.3.0

Template matching (`OrientationPlan`, `OrientationMatcher`), symmetry reduction, IPF colour keys and the comparison with pyxem
and orix are slice 5. A kinematic match cannot tell +c from −c on a polar structure: no consumer may claim a polarity from it.
