# DiffractoKit 0.0.3: what mac4DSTEM and macTEM change to import it

Written 2026-10-15 by DiffractoKit's session, for mac4DSTEM's and macTEM's sessions. `main` (56b5327 when this was
written, plus the consumer-check extension below) carries slices 1 to 3 and is what the tag `0.0.3` will be. The tag is
the owner's to make; until it exists nothing here can be pinned. The apps' sessions do the switch, each as a reviewed
tier-1 lane with its own registration wherever a number moves (the list is under "What numbers move"). The kit was read
against mac4DSTEM `7d6fdbd1` and macTEM `9482831` by `git show` only; recheck both files-and-lines at the bump.

What 0.0.3 holds, one product per module (all of them `public`, no `package` symbol is needed by an app):

| Product | Slice | Holds |
|---|---|---|
| `DiffractoCore` | 1 | `UnitCell` (all seven systems), `SymmetryOperation`, `IntegerMatrix`, `PointSymmetry`, `LaueClass`, `MillerIndex`, `Element`, `AtomSite`, `CrystalStructure`, `SymmetryDetection`, `CIFImport` / `CIFDocument`, `ReflectionList`, `Reflection`, `ReflectionFamily` |
| `DiffractoScattering` | 2 | `ElectronScatteringFactors` (Lobato and Van Dyck, py4DSTEM's table), `StructureFactors`, `StructureFactor`, `KinematicReflection`, `KinematicFamily`, `KinematicReflectionList` (`CrystalStructure.kinematicReflections`) |
| `DiffractoGeometry` | 3 | `ElectronWavelength`, `RingCalibration` (closed forms, `calibrate`, `FitReport`, `refusals`, `fit`), `RingFitCalibration`, `RingIndexing`, `ZoneAxis`, `DiffractionPattern`, `EllipseCalibration`, `EllipseCalibrationFit`, `EllipseProfileParameters`, `EllipseCalibrationModel`, `EllipseTransform` |
| `DiffractoOrientation` | later | a placeholder; do not link it yet |

`tools/consumer-check` (a separate package, the public API only, Debug, 0 warnings) passes 288 checks, 161 of them on
`DiffractoGeometry`. No negative controls yet (SpectroKit's `controls.py` pattern comes with the first tag).

## How to depend on it

- SwiftPM: `.package(url: "<the DiffractoKit remote, the owner's>", .exact("0.0.3"))`, then per target the products it
  imports: `.product(name: "DiffractoCore", package: "DiffractoKit")` and the same for `DiffractoScattering` and
  `DiffractoGeometry`. `DiffractoGeometry` depends on the other two, but `MemberImportVisibility` is on in both apps, so a
  file that names a type of a module imports that module itself.
- `mac4DSTEM.xcodeproj` needs the same product dependencies and build files for every target that links the code that
  moves, and its own pin (the `project.pbxproj` requirement and the workspace `Package.resolved`), as for SpectroKit.
- The package is macOS 27, Swift 6 language mode, tools 6.2, no unsafe flags: the app builds it optimised itself. The
  kit's tests run with `-c release` (`tools/test.sh`); an app's Debug build of a dependency is slow for the ring fit
  (`RingCalibration.fit` is a few hundred ms in release), so keep any `-O` workaround the apps already use for
  SpectroKit.
- Typed throws: `CIFImport.read` throws `CIFImportError`, `reflections(kMax:)` throws `ReflectionListError`,
  `kinematicReflections` throws `KinematicError`, `RingIndexing.index` throws `RingIndexing.Refusal`; a `catch` in an app
  that did `catch let e as X` still compiles.
- Names that exist in the apps and in the kit (`AtomSite`, `Reflection`, `DiffractionPattern`, `EllipseCalibration`,
  `EllipseCalibrationFit`, `RingCalibration`, `CIFImportError`) collide in a file that imports both: delete the app's
  copy in the same commit as the import, one type at a time.
- Licence: the kit is GPL-3.0 (`LICENSE`, `NOTICE`); both apps are GPL-3.0, so importing it is compatible. The kit's
  scattering table and the lattice conventions are py4DSTEM's (GPL-3.0); pymatgen (MIT) was the oracle and no code of it
  is copied.
- Pins behind the kit's numbers (for an app's own registration to cite): pymatgen 2026.9.24 (tag `v2026.9.24` =
  `8461938`) with pymatgen-core 2026.9.23 (`ae9a834`), spglib 2.8.0, NumPy 2.5.4, Python 3.14.7, SciPy 1.18.1 (the ring
  fit's `least_squares` optimum); py4DSTEM `f050d207851f8accd3aa550365baaf1f38230428` (`scattering_factors.txt`, sha256
  `be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4`, and the excitation-error formula); seeds
  mac4DSTEM `7d6fdbd1` and macTEM `9482831` (slice 3), macTEM `a58c84c` (slice 1).

## What numbers move

Read this before the file lists: these are the places where an app's output, a golden, a fixture or a stored string
differs after the switch. Each is a registered amendment of the kit, not an accident; each needs the app's own
registration (the old and the new value on the app's own fixture) in its tier-1 lane.

**Both apps (crystal, CIF, structure factors).**

1. **A CIF's symmetry loop no longer decides the Laue class.** The class is detected from the atoms (the loop's
   rotations are added to the candidates so a non-conventional setting is not under-detected; the loop is kept as
   `fileOperations`, and `fileOperationsAreSymmetries` says whether it was right). Both apps already confirmed
   cubic and hexagonal from the atoms, so a correct file reads the same; a file whose loop declares less than the atoms
   have (the apps' gold fixture lists one atom and the four centring translations) now reads its true class.
2. **The cell family is proposed at 1e-9 (lengths, relative) and 1e-9 degrees, not the apps' 1e-3 and 0.05 degrees**
   (`CIFImport.classifyFamily`). A cell the apps called cubic because its edges agree to 1e-3 (4.0782, 4.0785, 4.0782)
   is not cubic for the kit: a rotation must preserve every |g| to the ring rule's 1e-9, or an orbit splits across two
   rings. Such a CIF reads a lower Laue class (c/a = 1.0001 reads 4/mmm, 16 operations, as spglib), and the app's
   `ACOMCrystalSymmetry` for it becomes `.identity`. `UnitCell.metricFamily(lengthTolerance:angleToleranceDeg:)` takes
   the apps' tolerances if an app wants them; passing them is a decision to register, not a default.
3. **Trigonal and other lower-Laue structures import instead of being refused.** quartz (P3_221), calcite (R-3c),
   corundum, apatite (6/m), pyrite (m-3), a tetragonal structure with an accidentally cubic metric: the apps' importers
   throw `symmetryNotSupportedByStructure` (and `symmetryUnresolvableAtWrittenPrecision`); the kit reads -3 or -3m,
   6/m, m-3, 4/mmm and so on. The apps map their family from the class (below), and every class that is not m-3m,
   6/mmm or 2/m with b unique is `.identity` ("Unreduced"), which is what the apps already do for every other cell.
   This is a new accepted input, not a changed number; its UI wording (the apps' refusal texts) goes.
4. **`x/2`, `xy`, `x2`, `x.5` and singular rotations are refused.** The symmetry-operation grammar is mac4DSTEM
   `7d6fdbd1`'s `parseSymmetryComponent` / `parseFraction` (empty fraction pieces, a non-integer coefficient on x, y, z,
   a term that ends at a sign or the end), plus a refusal of a rotation with determinant other than +-1 (`x,x,z`), which
   neither app has. **macTEM at `9482831` still has the old reading**: its parser reads `x/2` as x + 2 (an empty
   numerator is dropped by `split`), so a file macTEM accepted with a wrong operation is now refused; mac4DSTEM's only
   new refusals are the singular rotation and, at the parser, a non-finite or over-1000 coefficient (its incomplete-list
   check guards that later).
5. **Magnesium (P6_3/mmc) has 25 symmetry-allowed rings before the |F| cut, not 24** (ideal hcp 24, not 23): {3 0 -3 1} is
   allowed by the operations and has F = 0 by the basis, exactly as Si's {2 2 2}. After the cut nothing moves: the apps'
   lists are always after the cut (Au 8, Al 8, Si 10, Mg 24, ideal hcp 23, beta'' 417 rings and 1442 reflections at
   1.2 A^-1), and `kinematicReflections(kMax:).families` reproduces those counts. A call that asks the kit for
   `ReflectionList.families()` (slice 1, before the cut) sees the 25th.
6. **Ring names.** Silicon's {5 1 1} + {3 3 3} ring is named `"{5 1 1} + {3 3 3}"`: the orbits of a ring run in
   descending lexicographic order of their canonical index (macTEM's `orbits` are `[(5 1 1), (3 3 3)]` too). Miller-Bravais
   (four-index) names appear only for a hexagonal-axes Laue class in a hexagonal metric (corundum in rhombohedral axes
   keeps three indices). A stored or displayed ring name, a golden string or a legend that compares names is checked
   against these.
7. **F is in A^-2, not "A".** F = (1/V) sum f_e e^(2 pi i h.r) with f_e in A and V in A^3, so F is A^-2, |F|^2 is
   A^-4 and the cut 1e-4 is A^-2. The apps' comments say "A" (py4DSTEM's name for f_e's unit). **The numbers are
   identical, bit for bit** (|g|, Re F, |F|^2 on the apps' presets; Im F below 1e-15 absolute is rounding noise of a
   centrosymmetric sum). Only comments, labels and any unit string change. The cut is a fixed number against a quantity
   that scales as 1/V, so a large cell keeps fewer reflections (a 30 A cell at 2 A^-1 keeps 16,326 of 904,088); that was
   always so and is now stated.
8. **Order of equal-|g| reflections is defined.** The list is sorted by |g| and reflections within 1e-12 relative of each
   other come in descending (h, k, l); the apps' `sort` on `gLength` alone was unstable between equal values. A consumer
   that is order-sensitive within a ring (the golden of a template, a first-found tie-break) may see a different but now
   reproducible order.
9. **Built through `CIFImport`, beta'' sums in another site order** than the apps' hand list: Re F differs by at most
   2.7e-13 and |F|^2 by 5.5e-13 of max(|F|, 1e-4) (measured; |g|, V and the kept sets are bit equal). Built from the
   apps' own hand site list the sums are bit identical at 1.2 and 2.0 A^-1. An app that moves its `betaDoublePrime`
   preset to `CrystalStructure` keeps the hand list to keep its goldens.
10. **A written 0.333 is not snapped to 1/3.** As in the apps (pymatgen does snap within 1e-4). With 1/3 written as 0.333,
    magnesium keeps 24 operation-absent reflections with |F| up to 3e-4 (306 kept and 28 rings against 270 and 24);
    `keptButOperationAbsent` names them. The apps keep them too (bit-for-bit kept sets), so nothing moves; the list can
    now say which ones they are.
11. **Smaller refusals.** A list is refused when the cell is invalid or the expected count (4/3) pi k_max^3 V exceeds ten
    million (`cellTooLargeOrInvalid`); a coplanar cell (volume under 1e-6 abc) is `degenerateVolume`; a NaN or infinite
    form factor, occupancy or coordinate is `invalidStructure`, never zero reflections; k_max above 4 A^-1 is refused.

**mac4DSTEM only (wavelength, ellipse).**

12. **The wavelength moves.** `DPC.electronWavelengthAngstrom(voltageKV:)` (12.2639 / sqrt(V (1 + 0.97845e-6 V))) is
    2.6e-5 to 2.8e-5 relative below the kit's CODATA-2018 lambda (80 kV: -2.84e-5; 120 kV: -2.80e-5; 200 kV: -2.72e-5;
    300 kV: -2.64e-5), and `ParallaxPreprocessor.electronWavelengthAngstrom(energyEV:)` (rounded constants) is 1.1e-7 to
    1.3e-7 below it (-1.23e-7 at 200 kV). The kit's lambda at 200 kV is 0.025079340450548003 A; the apps' are
    0.02507865865 (DPC) and 0.02507933736 (Parallax). What moves in mac4DSTEM by exactly that relative amount: the ACOM
    excitation error s_g = g.n + lambda |g|^2 / 2 (`AppState+ACOM` passes `DPC.electronWavelengthAngstrom` as
    `wavelengthAngstrom` to `OrientationPlan.generate`), so every curved template changes in its curvature term;
    `milliradiansPerDetectorPixel` (1000 lambda dq) in the calibration readout and the exported units; and every Parallax
    quantity that multiplies or divides by `calibration.wavelengthAngstrom` (probe angles in mrad, the aberration
    frequencies, the chi of depth sectioning and the aberration fit). A golden of any of these is re-recorded with the
    old and new values in the registration. The DPC form returns `nil` for a non-positive or non-finite voltage; the kit
    returns NaN, so a wrapper keeps the app's guard.
13. **`EllipseCalibration` and `Crystal.reflections` are bit identical to `7d6fdbd1`'s** on 35 rings x six calls (every
    output field, every error and fallback reason, at -O and at -Onone, each build held to the apps built the same way),
    and on the apps' presets and twelve random oblique cells: nothing moves here. `fitAmorphousRing`'s annulus is
    `inner < r < outer`, `fit1D`'s `inner < r <= outer`, kept for identity.

**macTEM only (ring fit, indexing).**

14. **The ring fit moves to mac4DSTEM `7d6fdbd1`'s (ADR 069).** macTEM's `EllipseCalibration` copy at `9482831` is
    the pre-ADR-069 file; `RingCalibration.fit` runs the kit's, with ADR 069's two starts and the 0.30 px ring-width
    floor. Clean plants move about 1e-9 relative; macTEM's S3-b noisy-ring numbers (the `RingFitCalibrationTests` bands
    and the harness numbers measured on its copy) will move and are re-measured and re-registered by macTEM; its eight
    `RingFitCalibrationTests` pass unchanged on the kit's fit. ADR 069's rates (2009 synthetic rows: wrong 414 -> 79,
    right 1401 -> 1814) are mac4DSTEM's, on its harness sets, not claimed for macTEM's.
15. **The closed forms do not move.** `RingCalibration.wavelengthAngstroms`, `twoTheta`, `cameraLengthPixels`,
    `qPerPixel`, `calibrate` and `refusals` are macTEM's, formula, constant, NaN rule, refusal order and message kept
    (the lambda is macTEM's CODATA-2018 one, so nothing moves for macTEM's wavelength). Au {111} at R 150 px, 200 kV:
    L_px 14082.307097295722, dq 0.0028313372405787253 A^-1/px, L 197.15229936214011 mm at 0.014 mm. Two refusals are
    additions already in macTEM's file (`braggLimit`, `invalidInput`); a non-finite fit number gives `[.invalidInput]`
    alone. One correction of a comment: the exact L = R / tan 2 theta is SMALLER than the small-angle R d / lambda by about
    1.5 (lambda/2d)^2 (4.25e-5 at 200 kV Au {111}), not larger.
16. **Indexing verdicts are reproduced exactly** (macTEM's T1 tables, T4 refusals and its 41 named cases, re-expressed on
    the kit's lists), and they now work for every crystal system: `unsupportedLattice` (macTEM indexed only cubic,
    hexagonal and b-unique monoclinic) disappears, a tetragonal, orthorhombic, trigonal or triclinic cell is indexed.
    `unsupportedElement` is no longer an indexing refusal: the list builder refuses it
    (`KinematicError.unsupportedElements`) before `index` is called. A free `kMax` is a parameter of the families
    overloads; the `KinematicReflectionList` overload uses the list's own, which is the one a caller should use
    (a larger k_max than the list's turns `beyondKMax` off).

## mac4DSTEM

Seed `7d6fdbd1`. Files it can drop or reduce:

| File | After the switch |
|---|---|
| `Core/Crystal/Crystal.swift` | Dropped except the presets and the app's cell constructors: `fcc/bcc/sc/diamond/hcp(a:z:)`, `Crystal.cubic(_:a:z:)`, `CubicStructure`, `aluminum ... magnesium`, `betaDoublePrime`, `tungstenDisulfide` become functions returning `CrystalStructure` (`CrystalStructure.detectingSymmetry(cell:sites:)`). `reflections(kMax:tolerance:)`, `latReal`, `latInv`, `volume`, the lattice algebra, `axisDelta`, `unsupportedElements` and `Reflection` go. |
| `Core/Crystal/CrystalModel.swift` | Only the app-side model stays: `id`, `displayName`, `source`, `symmetry: ACOMCrystalSymmetry`, `extraProvenance`, `provenance`, `contentFingerprint`, `revisionID`, `CrystalModelSelection`, `CrystalModelLibrary`, `customCubic`. The field `crystal: Crystal` becomes `structure: CrystalStructure`. `validationIssues` keeps the ACOM and model checks; the cell, occupancy, coordinate and close-contact checks are `CrystalStructure.validationIssues` and `shortestCloseContact` (check the issue codes and texts against the app's: they are the kit's). |
| `Core/Crystal/ScatteringFactors.swift` | Dropped: `ElectronScatteringFactors.f(z:gSquared:)` (table of 103 elements, identical rows), `isSupported(z:)`, `maximumAtomicNumber`. `ScatteringFactors.symbols` becomes `Element.symbol(forAtomicNumber:)` (Z up to 118). |
| `Core/Crystal/CIFImport.swift` | The parser, tokenizer, expansion, symmetry operations, `spaceGroupOrder`, `requiredOperationCount`, `classifyFamily`, `verifyFamily`, `confirmedFamily` and the position tolerances go: `CIFImport.read(_:fileBaseName:)` returns a `CIFDocument` (`structure`, `displayName`, `dataBlockName`, `declaredSpaceGroupSymbol`, `declaredSpaceGroupNumber`, `coarsestHalfStep`, `hadSymmetryOperations`). Stays: the model id (`imported_<sanitised base>`), `sanitize`, the mapping to `ACOMCrystalSymmetry` and `CrystalModel`. `displayName(dataBlockName:fileBaseName:)` is in the kit with the VESTA placeholder rule. |
| `Core/Analysis/EllipseCalibration.swift` | Dropped entirely (`EllipseCalibration`, `EllipseCalibrationFit`, `EllipseProfileParameters`, `EllipseCalibrationModel`, `FitError`, with the same names and the same labels). |
| `Core/Data/DiffractionPattern.swift` | `DiffractionPattern(qy:qx:pixels:)` stays as the display type (`FloatImage`, the contrast and the -1 sentinel are presentation); the fit takes the kit's three-field `DiffractoGeometry.DiffractionPattern(width: qx, height: qy, pixels:)`. The fit's `centerQX` is the ROW and `centerQY` the COLUMN, as now. |
| `Core/Data/Calibration.swift` | `ellipseTransform()`, `acceptedEllipse`, `ellipseCorrectedOffset`, `ellipseUncorrectedOffset` reduce to calls of `EllipseTransform(a:b:theta:)`, `.accepted`, `.correctedOffset(dx:dy:)`, `.uncorrectedOffset(dx:dy:)` (the same formulas, py4DSTEM `_transform`, `a < b` allowed). The stored `ellipseA/B/Theta`, `hasEllipse` and the provenance stay. |
| `Core/Analysis/DPC.swift`, `ParallaxPreprocessing.swift` | `electronWavelengthAngstrom(voltageKV:)` and `electronWavelengthAngstrom(energyEV:)` become one-line calls of `ElectronWavelength.angstroms(volts:)` (kV x 1000, eV as volts); the DPC form keeps its `nil` for a non-positive or non-finite voltage. **This is the number-moving change of item 12.** |
| `Core/Crystal/OrientationPlan.swift`, `PhaseReferenceLibrary.swift`, `Core/Analysis/OrientationResult.swift` (the geometry parts only) | `PhaseReferenceLibrary.lowIndexZoneAxes` -> `ZoneAxis.lowIndex` (the same 49 axes by value and in order); `cartesianZoneAxis(_:crystal:)` -> `ZoneAxis.cartesian(in: cell)`; `ACOMOrientation.detectorBasis(zoneAxis:)` -> `ZoneAxis.detectorBasis(_:)` (same formula); the excitation error and the (x, y) of `OrientationPlan.project` and `projectedVectors` -> `ZoneAxis.project(_:onto:wavelengthAngstroms:sgMax:)`. The weighting (intensity x Gaussian of s_g, `intensityPower`, `radialPower`, polar conversion) stays in the app. |

Type mapping:

| mac4DSTEM | Kit |
|---|---|
| `Crystal` (a, b, c, alpha, beta, gamma, `latReal`, `latInv`, `volume`) | `CrystalStructure` with `cell: UnitCell` (`UnitCell(a:b:c:alphaDeg:betaDeg:gammaDeg:)`, `realLattice`, `reciprocalLattice`, `volume`) |
| `AtomSite(z:fractional:occupancy:)` | `AtomSite(z:fractional:occupancy:label:)` |
| `Crystal.reflections(kMax:tolerance: 1e-4)` | `structure.kinematicReflections(kMax:cut: 1e-4).kept` (`StructureFactors.defaultCut`) |
| `Reflection` (`h, k, l`, `g`, `gLength`, `structureFactorRe`, `structureFactorIm`, `intensity`) | `KinematicReflection`: `reflection.hkl.h/.k/.l`, `reflection.g`, `reflection.gLength`, `structureFactor.real`, `.imaginary`, `.intensity` (the allowed flag and `d` are extra) |
| `Crystal.unsupportedElements` | `KinematicError.unsupportedElements([Int])`, thrown by `kinematicReflections` |
| `CIFImport.crystalModel(from:fileBaseName:)` | `CIFImport.read(_:fileBaseName:)`, then the app's model from `document.structure` |
| `CIFImportError` | `missingTag`, `unparsableNumber`, `noAtomSites`, `unknownElementSymbol`, `unsupportedSymmetryOperation`, `missingSymmetryOperations`, `truncatedLoop`, `multipleStructures`, `incompleteSymmetryOperations`, `symmetryOperationsOutsideStructure` keep their names and payloads; `invalidModel([String])` -> `invalidStructure([String])`; `unsupportedElement` -> `KinematicError.unsupportedElements`; `unsupportedPointGroup`, `symmetryNotSupportedByStructure`, `symmetryUnresolvableAtWrittenPrecision` have no counterpart (item 3) |
| `ACOMCrystalSymmetry` | **stays in the app.** The kit gives `structure.laueClass` (`LaueClass?`) and `laueGroup: PointSymmetry`; the app maps `.mThreeBarM` -> `.cubic`, `.sixOverMmm` -> `.hexagonal`, `.twoOverM` with `structure.cell.metricFamily() == .monoclinic(uniqueAxis: .b)` -> `.monoclinic`, everything else (and nil) -> `.identity`. The unique-axis test is needed: the kit's 2/m is any unique axis, ACOM's `.monoclinic` is b unique only. |
| `lowIndexZoneAxes: [SIMD3<Int>]` | `ZoneAxis.lowIndex: [ZoneAxis]` (`ZoneAxis(u, v, w)`, `.u/.v/.w`); a phase's own `zoneAxes` convert |
| `ACOMOrientation.detectorBasis(zoneAxis:)` | `ZoneAxis.detectorBasis(_:)` (the argument must be a unit vector: `ZoneAxis.isUnitDirection`) |
| `OrientationPlan.project`'s sg, x, y | `ZoneAxis.Projected` (`reflection`, `sg`, `x`, `y`) from `ZoneAxis.project`; it keeps nothing for a non-unit n, a non-finite or negative lambda |
| `DPC.electronWavelengthAngstrom(voltageKV:)` | `ElectronWavelength.angstroms(volts:)` |
| `EllipseCalibration.*` | the same names in `DiffractoGeometry` |
| `Calibration.ellipseTransform` / `ellipseCorrectedOffset` | `EllipseTransform` |
| `DiffractionPattern(qy:qx:pixels:)` | `DiffractionPattern(width: qx, height: qy, pixels:)` |

mac4DSTEM's tests: the app's `CrystalMathTests`, `CIFImportTests` and friends, `EllipseNoisyRingTests`,
`EllipseStopRuleTests` and `CalibrationEllipseRow29Tests` keep running unchanged against the kit-backed types as the
app's own check that the bit-identity claims hold in its build; the kit holds the same plants (SplitMix64 and Box-Muller
renderer, identical Float pixels) as its own gate.

Also in mac4DSTEM, not in 0.0.3: `PhaseVectorMatching`, the reference-library build, the orientation plan's templates and
the ACOM symmetry reduction (slice 4 and `DiffractoOrientation`), `QCalibration`, `CalibrationUnitConversion`, the
annulus guard (never landed). `MaterialsProjectImport` stays app-side: it calls `CIFImport.classifyFamily`,
`elementSymbol(from:)` and `elementSymbolToZ`, which go; the family classification there becomes the same Laue-class
mapping and the element lookups become `Element.normalizedSymbol(fromCIF:)` and `Element.atomicNumber(forSymbol:)`
(verify the strings it feeds them at the switch).

## macTEM

Seed `9482831`. Files it can drop or reduce:

| File | After the switch |
|---|---|
| `Crystal/Crystal.swift`, `Crystal/ScatteringFactors.swift`, `Crystal/CIFImport.swift` | As mac4DSTEM's: the copies go, the presets (`Crystal.fcc ... magnesium`, `betaDoublePrime`, `tungstenDisulfide`) become `CrystalStructure` builders. |
| `Crystal/CrystalModel.swift` | Only the app-side model stays (as mac4DSTEM's). |
| `Analysis/EllipseCalibration.swift` | Dropped (the pre-ADR-069 copy). |
| `Analysis/RingCalibration.swift`, `RingCalibration+Fit.swift` | Dropped: `RingCalibration`, `RingFitCalibration`, `RingCalibration.fit` are the kit's with the same names, labels and refusals. At `9482831` nothing in the `macTEM` app target or in `Sources/TEMCore` outside these files names them: the callers are the tests (`RingCalibrationTests`, `RingFitCalibrationTests`, `DSpacingIndexingTests`, `MonoclinicIndexingTests`) and the harnesses (`tools/ring-calibration-parity`, `tools/indexing-parity`, `tools/ellipse-parity`), so the switch is the imports plus the pattern conversion. |
| `Analysis/DSpacingIndexing.swift` | Dropped: `lattice(of:)`, `orbit(of:in:)`, `Family`, `FamilyList`, `families(of:kMax:)`, `HKL`, `Lattice` and `unsupportedLattice` go; `index` and its types are `RingIndexing`. |
| `Data/DiffractionPattern.swift` | Stays as the app's display type; the fit takes `DiffractoGeometry.DiffractionPattern(width: qx, height: qy, pixels:)`. |

Type mapping:

| macTEM | Kit |
|---|---|
| `Crystal`, `AtomSite`, `Reflection`, `ScatteringFactors`, `CIFImport` | as mac4DSTEM's |
| `DSpacingIndexing.Family` (g, `members: [HKL]`, `orbits: [HKL]`, intensity, lattice) | `KinematicFamily` (`family: ReflectionFamily` with `g`, `members: [MillerIndex]`, `orbits`, `name`, `d`, `multiplicity`; `intensity` = sum of |F|^2 over the members) |
| `DSpacingIndexing.FamilyList(lattice:kMax:families:)` from `families(of: crystal, kMax:)` | `KinematicReflectionList` from `structure.kinematicReflections(kMax:)` |
| `DSpacingIndexing.index(q:sigma:in: FamilyList, coverage:)` | `RingIndexing.index(q:sigma:in: KinematicReflectionList, coverage:)` (the list's own k_max); `index(q:sigma:in: [KinematicFamily], kMax:coverage:)` and `[ReflectionFamily]` overloads for a bare list |
| `DSpacingIndexing.Assignment`, `Candidate` (`best`, `next`, `margin`, `z`, `deviation`, `relativeDeviation`, `dMeasured`, `validation`) | `RingIndexing.Assignment`, `RingIndexing.Candidate` with the same fields (`family` is a `ReflectionFamily`) |
| `DSpacingIndexing.Refusal` (`invalidInput`, `unsupportedLattice`, `unsupportedElement`, `ambiguous`, `beyondKMax`, `noneWithinTolerance`) | `RingIndexing.Refusal` (`invalidInput`, `ambiguous`, `beyondKMax(reach:kMax:)`, `noneWithinTolerance(nearest:)`); the other two cannot occur (item 16) |
| `DSpacingIndexing.defaultCoverage`, `minimumMargin`, `reach(coverage:)`, `validation` | the same names on `RingIndexing` |
| `RingCalibration.*`, `RingFitCalibration` | the same names in `DiffractoGeometry` |
| `EllipseCalibration.*` | the same names in `DiffractoGeometry` (the ADR-069 file, item 14) |
| `ACOMCrystalSymmetry` (`Analysis/OrientationResult.swift`) | stays in the app; mapped from `LaueClass` as in mac4DSTEM's table |

macTEM's tests and tools: the harnesses compile the app's sources through `tools/lib/copies.manifest` (mac4DSTEM's through `tools/lib/sources.manifest`), so each manifest line for a dropped file is replaced by the kit's product at the switch, as SpectroKit's was. `tools/indexing-parity` and `tools/ring-calibration-parity` (the 320-case decimal reference,
68 refusal plants, 18 Bragg-limit plants, 10 radius plants) are the kit's pins too (ported to `tools/slice-3/`), so macTEM
can keep running them against the kit-backed types as its parity gate; `tools/crystal-parity` is slice 2's reference.

## Order of work the kit suggests

Not binding; the apps' sessions decide. Each step is its own commit and its own tier.

1. **mac4DSTEM, wavelength** (tier 1, small): the three-line replacement and the golden re-record of item 12 with the old
   and new numbers in the registration.
2. **Both, the crystal types** (tier 1 for the CIF/Laue/refusal changes, items 1 to 11): `Crystal` -> `CrystalStructure`,
   `CIFImport` -> `CIFImport.read`, `ScatteringFactors` -> `ElectronScatteringFactors`; the apps' own CIF fixtures
   (gold, silicon, magnesium, beta'', the seeded quartz, calcite and pyrite) are the registration's cases. The
   kit holds completed versions of the seeded quartz, calcite and pyrite (their operation loops are incomplete as seeded:
   the apps used them to show a refusal, not a cell content), marked `# MODIFIED`.
3. **mac4DSTEM, ellipse** (tier 1 but bit identical): `EllipseCalibration`, `DiffractionPattern` conversion,
   `EllipseTransform`.
4. **macTEM, ring calibration and indexing** (tier 1; item 14 moves its numbers).
5. **mac4DSTEM, zone axes** (tier 2 for `lowIndexZoneAxes`, `cartesian`, `detectorBasis`; tier 1 where the wavelength
   of step 1 reaches `project`).

The tag is the owner's. After it, `Package.swift` and the Xcode project pin `.exact("0.0.3")`; a later tag is a new note.

## Not mapped / known limits

- No kit counterpart: the app presets and the cubic-structure editor (`Crystal.CubicStructure`, `Crystal.cubic`), the model
  library, `MaterialsProjectImport` (the family classifier it borrows is gone), the ACOM symmetry reduction, fundamental
  zones and IPF keys (`DiffractoOrientation`, empty), `PhaseVectorMatching` and the reference library, `QCalibration`,
  `CalibrationUnitConversion`, the display helpers of `DiffractionPattern`, `FloatImage`, `RGBAImage`.
- Validation: none of it is validated on a real diffraction pattern. `RingIndexing`'s `accepted` means
  "19 times as likely as the next listed family at equal priors under a Gaussian sigma_q", not "the ring is this family"
  (`RingIndexing.validation` says so on every assignment); macTEM's T3 and T5 numbers stand. No real CIF from a database
  (COD, ICSD) has been read yet by the kit.
- `ZoneAxis.project` keeps py4DSTEM's overall sign deviation: n is the direction the beam TRAVELS, s_g = g.n + lambda
  |g|^2 / 2; py4DSTEM's is -s_g / (1 + lambda g.n) (both terms flip, the zero locus is the same) and the denominator is
  not applied. mac4DSTEM's own comment says only the curvature term changes sign; that wording is the kit's correction.
- A primitive (non-conventional) P1 cell of a centred lattice detects a lower Laue class than its conventional one (slice
  1b is the lattice reduction); with the file's operation loop the rotations join the candidates and these cells agree.
- Pseudo-symmetry within the 5e-3 fractional tolerance reads as symmetry (0.25 A in a 50 A cell); a Cartesian tolerance
  is a candidate for slice 1b.
