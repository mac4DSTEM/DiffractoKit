# Changelog

All notable changes to DiffractoKit. Versions are semantic (0.x: a minor may change the API). The numbers each module
reproduces, and against what, are in `docs/parity.md`.

## Unreleased (0.3.0)

Vector-matched phase mapping and the value types of an orientation result, in `DiffractoOrientation`.

### Added
- `PhaseReferenceLibrary` (`build(phases:settings:)`, `Failure`, `lowIndexZoneAxes` (49 directions), `inPlaneSteps`,
  `projectedVectors`, `rotate`, `reflections(of:kMaxInvAngstrom:)`, `unsupportedElements(of:)`), `PhaseDefinition`,
  `ReferenceVector`, `PhaseOrientationReference` (`chanceMatchFraction`), `LatticeVector`, `OrientationRelationship`,
  `PhaseReferenceSettings` (7 settings).
- `PhaseVectorMatcher` (`map`, `classify`, `classifyKnownVariants`, `fitMatrixOrientation`, `fitZoneAxis`, `claims`,
  `experimentalVectors`, `nearest`, `score`, `containsFriedelPair`, `ZoneAxisFit`), `PhaseVectorSettings` (17 settings,
  `ClassificationRule`), `PhaseVectorResolution`, `PhaseVerdict`, `PhaseVectorResult`, `PhaseMap`, `PhaseDiskClaim`,
  `PeakCoordinates` (detector pixels with the origin subtracted in `Float`, or calibrated Å⁻¹).
- `EulerAngles` (Bunge z-x-z as SciPy, from and to py4DSTEM's orientation matrix), `ZoneAxisFrame` (`detectorBasis(zoneAxis:)`,
  `matrix(detectorBasis:inPlaneAngle:)`), `OrientationResult`.
- Every setting is a parameter with a documented default (the values tuned on the suite's and Thronsen et al.'s data); scans
  are `[[SIMD2<Double>]]`; progress and cancellation are closures; `map(maximumWorkers:)`.
- Regression tests: py4DSTEM f050d207 (`generate_diffraction_pattern` at λ → 0: hkl sets, |q| within 1.1e-15, the in-plane
  frame rotation per zone axis, intensities within py4DSTEM's complex64 precision), SciPy 1.18.1 Euler angles on 1000 rotations
  and the gimbal cases, NumPy closed forms, and the original Swift implementation bit for bit (about 1.68 million values).

### Fixed (against the implementation it was ported from)
- A rotation about the beam whose R22 rounds to 1 − 2⁻⁵³ (with R20 = R21 = 0) was lost by the Euler decomposition (both
  angles 0); it is folded into φ1.
- An exact score tie between two phases was broken by the process's hash seed; it goes to the lower phase index.

### Migration
For code moving from the app implementation: `Crystal` becomes `CrystalStructure` (`PhaseDefinition.crystal`); a structure the
kit refuses throws `PhaseReferenceLibrary.Failure.structureRefused`; `BraggVectors` becomes `scan:scanWidth:scanHeight:` with
`PeakCoordinates` (pass each peak's `Float` x, y as `SIMD2<Double>`; the pixel form reproduces the old arithmetic bit for bit);
the cancellation token becomes `progress:` / `isCancelled:` closures; `ACOMOrientation` is `ZoneAxisFrame`; `fitZoneAxis`
throws; a scan whose position count is not width × height gives nil (or `[]`) from every entry point, and an out-of-range
`matrixEntryIndex` makes `map` return nil. No other number moves.

## 0.2.0

The generality pass: every policy a parameter with the previous value as its documented default, the platform floor lowered,
the docs standing alone. The default path is bit identical to 0.1.1 (31,725 outputs compared).

### Added
- Policy types, each `Sendable, Equatable` with a defaulted memberwise `init` and `static let default`:
  `ReflectionList.Policy` (k_max cap 4 Å⁻¹, reflection and enumeration caps 1e7 and 1e8, the ring rule 1e-9, absence 1e-9),
  `CrystalStructure.ValidationPolicy` (close contact 0.5 Å, split-site occupancy 0.02), `SymmetryDetection.Policy`,
  `SymmetryOperation(xyz:maximumMagnitude:)` (1000), `CIFImport.Policy`, `EllipseCalibration.Policy` (18 fields plus
  `evidence`), `RingEvidence.Policy` (8), `RingCalibration.Policy`, `RingIndexing` `minimumMargin:` (ln 19). The |F| cut
  (`cut:`, 1e-4 Å⁻²) is documented as a policy. Where each is passed: the doc comment of each type.
- `kinematicReflections(..., retainingAll: Bool = true)`: with `false` the list from before the |F| cut is never held (`all` is
  empty; `kept` and `families` identical). A 30 Å cell at k_max 2.5 (1,767,026 reflections, 441,710 kept): 820 MB peak and
  2.3 s by default, 279 MB and 2.5 s with `false`.
- `CrystalStructure.symmetryRefusal`, `SymmetryDetection.Result.refusal`, `SymmetryDetection.Tolerance.refusal(for:)` and
  `isValid(for:)`, `ConventionalCell.reductionConverged`.
- A standalone README with the first ten minutes, `Examples/diffracto-example`, a DocC catalogue per module, CI.

### Changed
- The platform floor is macOS 15 and iOS 18 (was macOS 27); numbers unchanged there.
- No app or suite name remains in sources, test names or comments; test names changed accordingly.
- `ReflectionListError.cellTooLargeOrInvalid`'s text names the estimated count, the default limit and the policy field.

### Fixed
- Traps: `EllipseCalibration.Policy(angularBinCount:)` with a huge count; a symmetry position tolerance that is negative,
  denormal or tiny (now refused: the window must be 1e-9 up to 0.5 of a cell axis); a CIF with an absurd cell (overflow in the
  conventional-cell route: refused by validation or flagged `reductionConverged == false`); `IntegerMatrix.determinant` wraps.
- A position tolerance of 0 or NaN returned a structure without the identity operation; it is refused, and the identity is
  always present.
- A cell sheared by hundreds of units was silently left unreduced (slice 1b's review); a non-converged reduction is flagged.
- NaN, negative or infinite metric tolerances and Niggli epsilons are read as documented (0, or unreduced and unconverged).
- `RingIndexing.reach` answers NaN where `index` refuses the margin.

### Migration
Nothing is required: every new parameter defaults to the previous value. Notes for a caller who sets policies:
`EllipseCalibration.Policy(angularBinCount:)` derives the degeneracy and sparse bounds, while setting `angularBinCount` on
`.default` keeps 30 and 12; through `RingCalibration.fitRing` the calibration policy's evidence replaces the ellipse policy's,
and loosening coverage needs both `RingCalibration.Policy.minimumOccupiedAngularBins` and the ellipse policy's
`degeneracyBoundBins`; removing the `ReflectionList.Policy` caps removes the memory and time protection (a k_max of 100 reached
10 GB). Error texts quote the default bounds. A kinematic pattern cannot tell +c from −c on a polar structure: no consumer may
claim a polarity from a kinematic match. The app-specific notes of earlier releases are in git history
(`docs/consumers-0.0.3.md` to `docs/consumers-0.3.0.md`; they leave the repository at the next lean pass and stay in git history).

## 0.1.1

Slice 1b: Niggli reduction, the conventional cell in both settings, the Cartesian symmetry tolerance, detection on a cell list;
the CI type-check fixes.

## 0.1.0

Slices 1 to 3: cells of the seven systems, symmetry from the atoms, CIF import, reflections, absences, orbits and rings;
electron scattering factors and structure factors; the wavelength, ring and ellipse calibration, camera length, indexing and
zone axes; error texts; progress and cancellation of the fits.
