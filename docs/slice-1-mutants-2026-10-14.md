# Slice 1 break-once proof (2026-10-14 night)

Registration section 6, gate 3: every new test broken once by a mutant aimed at its rule. `tools/slice-1/mutants.py` applies each
planned mutant to a copy of the committed tree (`git archive HEAD`), runs `swift test --filter DiffractoCoreTests` (debug) and
requires red. Two runs, verbatim from the script's output (the scratch path shortened):

## Run 1, commit b4cfe17 (the frozen slice-1 commit): 16 of 19 red

The three green ones were not missing rules: floor vs ceil on an integer index bound is an equivalent mutant; loosening the
expansion's fractional tolerance alone is masked by its second guard (0.10 Å in real space); a ±0.002 dumbbell about a face
centre is cubic within the 5e-3 tolerance (every operation maps each half onto its own neighbourhood), so the
target-distinct guard never engaged. Each was replaced: the tile taken from one axis for all three, both merge guards
loosened, and the apps' bijection fixture (two extra atoms inside the tolerance of the origin atom).

```
mutants: 19 on commit b4cfe17, scratch <scratch>
baseline: green
RED                  reciprocal as G A, not G^-1 A [UnitCell.swift; aimed at UnitCellTests.testReciprocalRowsAreDualToTheRealOnes; caught elsewhere; 4 s]
RED                  beta read from the alpha tag [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 8 s]
GREEN (test missing) index tile floored instead of ceiled [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; not caught; 6 s]
RED                  orbit by the column action R hkl, not the row action hkl R [MillerIndex.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 17 s]
RED                  canonical key without the non-negative count [PointSymmetry.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 13 s]
RED                  inversion dropped from the monoclinic holohedry [PointSymmetry.swift; aimed at PointSymmetryTests.testHolohedriesAreTheRegisteredGroups; caught by it; 7 s]
RED                  C2 about a instead of b [PointSymmetry.swift; aimed at PointSymmetryTests.testSpecialOrbitSizes; caught by it; 6 s]
RED                  absence test without the R-fixes-h guard [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught elsewhere; 9 s]
RED                  family tolerance 1e-4 [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 9 s]
RED                  family tolerance 0 [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 10 s]
RED                  multiplicity as the orbit count [ReflectionList.swift; aimed at CIFImportTests.testGoldP1DetectsFm3mAndHasEightFamilies; caught by it; 9 s]
GREEN (test missing) expansion merge tolerance 0.5 [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; not caught; 20 s]
RED                  centring translation skipped in the expansion [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught by it; 15 s]
RED                  1/2+x read as x [SymmetryOperation.swift; aimed at SymmetryOperationTests.testScrewAxisParsesRotationAndTranslation; caught by it; 6 s]
RED                  Laue table: 8-with-4-fold branch inverted [PointSymmetry.swift; aimed at PointSymmetryTests.testLaueClassesOfTheSmallerGroups; caught by it; 6 s]
GREEN (test missing) detection without the target-distinct guard [SymmetryDetection.swift; aimed at SymmetryDetectionTests.testSplitSiteGoldIsNotCubic; not caught; 6 s]
RED                  detected translations not snapped to k/24 [SymmetryDetection.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 6 s]
RED                  the file's loop taken as the operations, not detected [CIFImport.swift; aimed at CIFImportTests.testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m; caught by it; 7 s]
RED                  proper-rotation orders read over every element [PointSymmetry.swift; aimed at CIFImportTests.testTrigonalInHexagonalAxesIsThreeBarM; caught by it; 7 s]
mutants: 16 of 19 red on commit b4cfe17
```

## Run 2, commit 9981f43: 19 of 19 red

```
mutants: 19 on commit 9981f43, scratch <scratch>
baseline: green
RED                  reciprocal as G A, not G^-1 A [UnitCell.swift; aimed at UnitCellTests.testReciprocalRowsAreDualToTheRealOnes; caught elsewhere; 8 s]
RED                  beta read from the alpha tag [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 10 s]
RED                  index tile from the second axis for every axis [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 5 s]
RED                  orbit by the column action R hkl, not the row action hkl R [MillerIndex.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 7 s]
RED                  canonical key without the non-negative count [PointSymmetry.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 5 s]
RED                  inversion dropped from the monoclinic holohedry [PointSymmetry.swift; aimed at PointSymmetryTests.testHolohedriesAreTheRegisteredGroups; caught by it; 4 s]
RED                  C2 about a instead of b [PointSymmetry.swift; aimed at PointSymmetryTests.testSpecialOrbitSizes; caught by it; 5 s]
RED                  absence test without the R-fixes-h guard [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught elsewhere; 10 s]
RED                  family tolerance 1e-4 [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 15 s]
RED                  family tolerance 0 [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 16 s]
RED                  multiplicity as the orbit count [ReflectionList.swift; aimed at CIFImportTests.testGoldP1DetectsFm3mAndHasEightFamilies; caught by it; 16 s]
RED                  expansion merge: both guards loosened (tolerance + 0.5, 100 Å) [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught by it; 21 s]
RED                  centring translation skipped in the expansion [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught by it; 16 s]
RED                  1/2+x read as x [SymmetryOperation.swift; aimed at SymmetryOperationTests.testScrewAxisParsesRotationAndTranslation; caught by it; 9 s]
RED                  Laue table: 8-with-4-fold branch inverted [PointSymmetry.swift; aimed at PointSymmetryTests.testLaueClassesOfTheSmallerGroups; caught by it; 7 s]
RED                  detection without the target-distinct guard [SymmetryDetection.swift; aimed at SymmetryDetectionTests.testSplitSiteGoldIsNotCubic; caught by it; 30 s]
RED                  detected translations not snapped to k/24 [SymmetryDetection.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 6 s]
RED                  the file's loop taken as the operations, not detected [CIFImport.swift; aimed at CIFImportTests.testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m; caught by it; 5 s]
RED                  proper-rotation orders read over every element [PointSymmetry.swift; aimed at CIFImportTests.testTrigonalInHexagonalAxesIsThreeBarM; caught by it; 5 s]
mutants: 19 of 19 red on commit 9981f43
```

## Runs 3 to 5, after the Gate B fixes (2026-10-15)

- **f629722, 28 mutants: 23 red.** Two greens were stale debug builds (a constant folded into another file's default
  argument is not recompiled when only its own file changes; `mutants.py` now touches every DiffractoCore source after
  each edit, and the metric mutant is applied at its use site), one was an incomplete mutant (the size guard has two
  clauses), and three rules had no test of their own (now: exact k/24 translations for 4-decimal coordinates, every
  detected operation a symmetry of the sites on every fixture, a shifted rock salt written to 5 decimals).
- **bb184c7, 29 mutants: 27 red.** The absence-tolerance mutant stayed green on every case: the least-squares
  translation's error ē − R ē has h · (ē − R ē) = 0 whenever h R = h, so the tolerance never fired. It was removed and its
  mutant replaced by "refinement dropped".
- **16ce14d, 29 mutants: 29 red.**

```
mutants: 29 on commit 16ce14d, scratch <scratch>
baseline: green
RED                  reciprocal as G A, not G^-1 A [UnitCell.swift; aimed at UnitCellTests.testReciprocalRowsAreDualToTheRealOnes; caught elsewhere; 4 s]
RED                  beta read from the alpha tag [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 10 s]
RED                  index tile from the second axis for every axis [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 8 s]
RED                  orbit by the column action R hkl, not the row action hkl R [MillerIndex.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 10 s]
RED                  canonical key without the non-negative count [PointSymmetry.swift; aimed at PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps; caught by it; 9 s]
RED                  inversion dropped from the monoclinic holohedry [PointSymmetry.swift; aimed at PointSymmetryTests.testHolohedriesAreTheRegisteredGroups; caught by it; 8 s]
RED                  C2 about a instead of b [PointSymmetry.swift; aimed at PointSymmetryTests.testSpecialOrbitSizes; caught by it; 8 s]
RED                  absence test without the R-fixes-h guard [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught elsewhere; 12 s]
RED                  family tolerance 1e-4 [ReflectionList.swift; aimed at PinsTests.testEveryStructureAgreesWithPymatgen; caught by it; 14 s]
RED                  family tolerance 0 [ReflectionList.swift; aimed at CIFImportTests.testBetaDoublePrimePredictedListAndFamilies; caught by it; 11 s]
RED                  multiplicity as the orbit count [ReflectionList.swift; aimed at CIFImportTests.testGoldP1DetectsFm3mAndHasEightFamilies; caught by it; 8 s]
RED                  expansion merge: both guards loosened (tolerance + 0.5, 100 Å) [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught by it; 9 s]
RED                  centring translation skipped in the expansion [CIFImport.swift; aimed at CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM; caught by it; 14 s]
RED                  1/2+x read as x [SymmetryOperation.swift; aimed at SymmetryOperationTests.testScrewAxisParsesRotationAndTranslation; caught by it; 11 s]
RED                  Laue table: 8-with-4-fold branch inverted [PointSymmetry.swift; aimed at PointSymmetryTests.testLaueClassesOfTheSmallerGroups; caught by it; 12 s]
RED                  detection without the target-distinct guard [SymmetryDetection.swift; aimed at SymmetryDetectionTests.testSplitSiteGoldIsNotCubic; caught by it; 11 s]
RED                  detected translations never put on the k/24 grid [SymmetryDetection.swift; aimed at RefuterRegressionTests.testCoarseCoordinatesGiveExactGridTranslations; caught by it; 9 s]
RED                  translation not refined (the trial translation kept) [SymmetryDetection.swift; aimed at RefuterRegressionTests.testShiftedRockSaltAtFiveDecimalsKeepsItsAbsences; caught elsewhere; 9 s]
RED                  the old snap: k/24 within the position tolerance [SymmetryDetection.swift; aimed at RefuterRegressionTests.testOneAtomAtAGeneralPositionIsPmThreeBarM; caught by it; 9 s]
RED                  final translation not re-verified [SymmetryDetection.swift; aimed at RefuterRegressionTests.testEveryDetectedOperationIsASymmetryOfTheSites; caught elsewhere; 9 s]
RED                  no closure: the operations returned as found [SymmetryDetection.swift; aimed at SymmetryDetectionTests.testSplitSiteGoldIsNotCubic; caught by it; 7 s]
RED                  images of any row merge [CIFImport.swift; aimed at RefuterRegressionTests.testCoSitedSpeciesAreNotMerged; caught by it; 10 s]
RED                  Miller-Bravais by the Laue class alone [ReflectionList.swift; aimed at RefuterRegressionTests.testRhombohedralAxesKeepThreeIndexNames; caught by it; 9 s]
RED                  no size guard [ReflectionList.swift; aimed at RefuterRegressionTests.testHugeCellsAreRefusedNotTrapped; caught elsewhere; 7 s]
RED                  a term not required to end at a sign (xy read as x+y) [SymmetryOperation.swift; aimed at SymmetryOperationTests.testRefusals; caught by it; 9 s]
RED                  singular rotations accepted [SymmetryOperation.swift; aimed at RefuterRegressionTests.testNonInvertibleOperationsAreRefused; caught by it; 9 s]
RED                  detection metric at the apps' 1e-3 [CrystalStructure.swift; aimed at RefuterRegressionTests.testNearlyCubicCellIsTetragonalAndNoRingIsNamedTwice; caught by it; 9 s]
RED                  the file's loop taken as the operations, not detected [CIFImport.swift; aimed at CIFImportTests.testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m; caught by it; 10 s]
RED                  proper-rotation orders read over every element [PointSymmetry.swift; aimed at CIFImportTests.testTrigonalInHexagonalAxesIsThreeBarM; caught by it; 8 s]
mutants: 29 of 29 red on commit 16ce14d
mutants exit 0
```
