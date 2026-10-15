# DiffractoKit 0.2.0: what mac4DSTEM and macTEM change to re-pin

Written 2026-10-15 by DiffractoKit's session, for mac4DSTEM's and macTEM's sessions. 0.2.0 is the generality pass: every
suite policy is a parameter with the old value as the documented default, no app or suite name stays in Sources, the
README and docs stand alone, and the platform floor is macOS 15 and iOS 18. Registered in `docs/registration-0.2.0.md`
(with its dated amendments). **TAG COMMIT: to be filled.** The tag is the owner's; until it exists nothing here can be
pinned. Read `docs/consumers-0.1.0.md` and `docs/consumers-0.1.1.md` first: this note lists only what changed since 0.1.1.

## In one table

| Who | Deletes or changes | Why |
|---|---|---|
| both | nothing is required | every new parameter defaults to the value the code used before; the default path is bit identical |
| both | the platform floor may be lowered to macOS 15 / iOS 18 | the kit imports only Foundation and simd; no availability branch; numbers unchanged at the floor |
| both | a decision per policy the app wants to own | a `Policy` value passed where listed below, instead of a fork or a wrapper |
| macTEM | its interim ring-fit bounds, if it held any beyond `RingCalibration.Policy` | `EllipseCalibration.Policy` and `RingCalibration.Policy` carry them |
| both | nothing for renamed tests | the kit's test names changed (no app or suite name in them); no app depends on them |

## What numbers move

None. The default path is bit identical: the full suite, the slice 1-3 pins, the 35 parity rings at -O and -Onone and the apps'
parity fixtures pass unchanged with the defaults. `Policy.default` equals the constants it replaces (a test per Policy type
checks every field against the old constant, and the mutants of `tools/kit-0.2.0/mutants-policy.py` change a default per
group). Error texts are unchanged: they still quote the default bounds (see "Not a parameter yet").

## The policies (additive)

Each type is `Sendable, Equatable`, has a memberwise `init` with every field defaulted, and `static let default`. A bound that is
NaN or negative is taken as documented on the type (refused or read as 0), never a trap. Where it is passed:

| Type | Fields (default) | Passed to |
|---|---|---|
| `ReflectionList.Policy` | `maximumKMax` 4 Å⁻¹, `maximumReflections` 1e7, `maximumEnumeration` 1e8, `familyTolerance` 1e-9 (the ring rule), `absenceTolerance` 1e-9 | `CrystalStructure.reflections(kMax:policy:)`, `ReflectionList.init(..., policy:)`, `kinematicReflections(..., policy:)`; the list keeps its `policy` |
| `CrystalStructure.ValidationPolicy` | `closeContactLimit` 0.5 Å, `splitSiteOccupancyTolerance` 0.02 (sum bound 1.02) | `validationIssues(using:)`, `isUsable(under:)`, `shortestCloseContact(using:)`, `kinematicReflections(..., validation:)` |
| `SymmetryDetection.Policy` | `occupancyTolerance` 0.02, `translationGrid` 1/24, `gridSnapFactor` 3, `precisionLimitFactor` 10 | `detectingSymmetry(..., policy:)`, `conventionalSetting`, `SymmetryDetection.detect(..., policy:)` |
| `SymmetryOperation(xyz:maximumMagnitude:)` | 1000 (`defaultMaximumMagnitude`): a larger number is unparsable | the initialiser |
| `CIFImport.Policy` | `maximumMergeDistance` 0.10 Å, `minimumMergeTolerance` 1e-6, `maximumCoordinateHalfStep` 1e-2, `rightAngleToleranceDeg` 0.05, `maximumOperationMagnitude` 1000, `positionTolerance` `.default`, `metricLengthTolerance`, `metricAngleToleranceDeg` (slice 1's), and the nested `symmetry` and `validation` policies | `CIFImport.read(_:fileBaseName:policy:)` |
| `StructureFactors` cut | 1e-4 Å⁻², already `cut:` since 0.0.1, now documented as a policy | `kinematicReflections(kMax:cut:...)` |
| `EllipseCalibration.Policy` | 18 fields: `angularBinCount` 36, `degeneracyBoundBins` 30, `sparseFloorBins` 12, `strongSignalFraction` 0.2, `oneRingRatio` 1.10, `minimumConicSamples` 20, `minimumProfileSamples` 32, `minimumRelativeSignal` 1e-6, `impulseRangeFactor` 2, `maximumConicResidual` 0.2, `maximumProfileResidual` 0.22, `profileSeedResidual` 0.4, `maximumAxisRatio` 5, `robustCutMultiplier` 4 · 1.4826, `minimumInlierPercent` 70, `ringWidthFloor` 0.25 px, `ringWidthFloorTolerance` 0.05, `loweredRingWidthFloor` 0.01 px; plus `evidence` | `policy:` on every ellipse and ring-profile fit; an invalid policy (`isValid` false) throws `FitError.invalidInput` |
| `RingEvidence.Policy` | `minimumSectorPixels` 8, `minimumAmplitudeFraction` 0.5, `minimumSignificance` 3, `minimumOccupiedBins` 30, `minimumDrop` 0.25, `maximumShift` 0.25 px, `shiftSignificance` 4, `maximumSecondaryPeak` 0.25 | `EllipseCalibration.ringEvidence(..., policy:)`, `isPeak(_:)`, `isSingleRing(_:)`, `geometryHolds(_:)`, `passes(_:)` |
| `RingCalibration.Policy` | `angularBinCount` 36, `minimumOccupiedAngularBins` 30, `maximumResidual` 0.22, `maximumEllipticity` 0.093, `evidence` | `RingCalibration.refusals(_:highTensionVolts:policy:)`, `fitRing(..., policy:ellipsePolicy:)`, `Refusal.reason(policy:)`; the fit uses the calibration policy's evidence |
| `RingIndexing` margin | `minimumMargin` ln 19 (`RingIndexing.minimumMargin`) | `index(...)` and `reach(coverage:minimumMargin:)` gained `minimumMargin:`; `coverage` was already a parameter |

An app that wants a stricter or looser rule sets one field on `.default` and passes the value; an app that passes nothing
changes nothing.

## Fixes since 0.1.1 (slice 1b's Gate B)

- **F1, multi-unit shears.** A cell sheared by hundreds of units, which 0.1.1 left silently unreduced, is now reduced (bit
  identical on the 2827 pinned cells), and a reduction that does not converge is flagged by `ConventionalCell.reductionConverged` (and
  `ConventionalSetting.reductionConverged`) rather than passed on as reduced. Real files do not reach this.
- **F2-F4 are documented, not changed:** the monoclinic a <-> c relabelling, the angle tolerance not acting on the conventional
  route (pass `metricLengthTolerance` of 1e-5 or more for rounded primitive cells), `.cartesian(t)` acting like about t/2 of
  spglib's symprec. Read them in `docs/consumers-0.1.1.md`.

## Docs, names, CI

- The README stands alone (a CIF in; d-spacings and a ring calibration out; the README's program is the first example).
  `tools/readme-consumer` compiles only the README's code, in CI. A DocC catalogue per module with one article (`docc
  convert`, 0 warnings). The only name in the README is the repository's URL.
- No app or suite name remains in Sources, test comments or names, Package.swift or the consumer check's sources. Renamed tests
  do not affect an app that pins the kit. What keeps the names: the registrations, these notes, status, parity, NOTICE, the
  CIF fixtures' sha-pinned provenance lines.
- CI type-check: 25 expressions at or above 40 ms were rewritten into typed locals (worst 372 ms before); none is left above
  40 ms; values unchanged.

## Not a parameter yet

Constants that are policy-like but stayed constants in 0.2.0 (a later minor, or on an app's asking): the profile-fit
plausibility bounds, the cell-volume floor 1e-6, the Niggli epsilon (`UnitCell.niggliRelativeEpsilon`, already an argument of
`niggliReduced(relativeEpsilon:)`). Error texts quote the default bounds even when the caller passed another policy, except
`RingCalibration.Refusal.reason(policy:)`. The consumer check has no Policy checks yet. Not policies, and staying constants:
CODATA constants, the scattering table, closed forms, the Bragg-limit bound 1/sqrt(2), and numerical-method constants that
change no accepted result (the classification is in `docs/registration-0.2.0.md` §4).

## The consumer check

537/537 at the macOS 15 floor (Debug, 0 warnings), unchanged from 0.1.1; the iOS 18 build succeeds. The package suite is 306/306 at
6c16c9a before the policy tests (the policy tests bring 348 on their base; the merged count is filled at the freeze).

## mac4DSTEM

Re-pins with no code change. It may now hold its own ring-fit and symmetry choices as `Policy` values (its lane's decision),
and may lower its floor to macOS 15 if it wishes; the app stays at 27.

## macTEM

Re-pins with no code change. Its Diffraction room can expose `RingCalibration.Policy` bounds without a fork. The tag is the
owner's; re-run the ring-fit parity and the real-pattern measurements after the re-pin, as for 0.1.0.
