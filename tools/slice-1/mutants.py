#!/usr/bin/env python3
"""Break-once proof for slice 1 (registration section 6, gate 3): each planned mutant is applied to a COPY of the
committed tree (git archive HEAD into a scratch directory; the repository is never edited), DiffractoCoreTests is run,
and the run must go RED. A mutant that stays green means a test is missing. Prints one line per mutant and a summary;
exits non-zero unless every mutant is red.

    python3 -I tools/slice-1/mutants.py [<scratch-dir>]        (default: <repo>/.build/mutants)

Each mutant is (name, file, exact old text, new text, the test it is aimed at); the old text must occur exactly once.
"""
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRATCH = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, ".build", "mutants")
SRC = "Sources/DiffractoCore/"

MUTANTS = [
    ("reciprocal as G A, not G^-1 A", SRC + "UnitCell.swift",
     "self.reciprocalLattice = UnitCell.matMul(metricInverse, real)",
     "self.reciprocalLattice = UnitCell.matMul(metric, real)",
     "UnitCellTests.testReciprocalRowsAreDualToTheRealOnes"),
    ("beta read from the alpha tag", SRC + "CIFImport.swift",
     'let beta = try requiredNumber("_cell_angle_beta")',
     'let beta = try requiredNumber("_cell_angle_alpha")',
     "CIFImportTests.testBetaDoublePrimePredictedListAndFamilies"),
    # Floor instead of ceil is an equivalent mutant (|h| ≤ k_max |a| for an integer h is |h| ≤ floor(k_max |a|)); the
    # real failure mode is a bound taken from one axis for all three, which loses reflections in an oblique cell.
    ("index tile from the second axis for every axis", SRC + "ReflectionList.swift",
     "let tile = cell.realLattice.map { max(0, Int((kMax * length($0)).rounded(.up))) }",
     "let tile = cell.realLattice.map { _ in max(0, Int((kMax * length(cell.realLattice[1])).rounded(.up))) }",
     "PinsTests.testEveryStructureAgreesWithPymatgen"),
    ("orbit by the column action R hkl, not the row action hkl R", SRC + "MillerIndex.swift",
     "h * rotation.rows[0].x + k * rotation.rows[1].x + l * rotation.rows[2].x,",
     "h * rotation.rows[0].x + k * rotation.rows[0].y + l * rotation.rows[0].z,",
     "PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps"),
    ("canonical key without the non-negative count", SRC + "PointSymmetry.swift",
     "(lhs.nonNegativeCount, lhs.h, lhs.k, lhs.l) < (rhs.nonNegativeCount, rhs.h, rhs.k, rhs.l)",
     "(0, lhs.h, lhs.k, lhs.l) < (0, rhs.h, rhs.k, rhs.l)",
     "PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps"),
    ("inversion dropped from the monoclinic holohedry", SRC + "PointSymmetry.swift",
     "case .b: generators = [c2y, .inversion]",
     "case .b: generators = [c2y]",
     "PointSymmetryTests.testHolohedriesAreTheRegisteredGroups"),
    ("C2 about a instead of b", SRC + "PointSymmetry.swift",
     "case .b: generators = [c2y, .inversion]",
     "case .b: generators = [c2x, .inversion]",
     "PointSymmetryTests.testSpecialOrbitSizes"),
    ("absence test without the R-fixes-h guard", SRC + "ReflectionList.swift",
     "for op in operations where hkl.applying(op.rotation) == hkl {",
     "for op in operations {",
     "CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM"),
    ("family tolerance 1e-4", SRC + "ReflectionList.swift",
     "public static let familyTolerance: Double = 1e-9",
     "public static let familyTolerance: Double = 1e-4",
     "PinsTests.testEveryStructureAgreesWithPymatgen"),
    ("family tolerance 0", SRC + "ReflectionList.swift",
     "public static let familyTolerance: Double = 1e-9",
     "public static let familyTolerance: Double = 0",
     "CIFImportTests.testBetaDoublePrimePredictedListAndFamilies"),
    ("multiplicity as the orbit count", SRC + "ReflectionList.swift",
     "public var multiplicity: Int { members.count }",
     "public var multiplicity: Int { orbits.count }",
     "CIFImportTests.testGoldP1DetectsFm3mAndHasEightFamilies"),
    # The merge has two guards (fractional tolerance AND 0.10 Å in real space); loosening one alone is masked by the
    # other, so the mutant loosens both: distinct sites then merge and the cell content shrinks.
    ("expansion merge: both guards loosened (tolerance + 0.5, 100 Å)", SRC + "CIFImport.swift",
     ["let tolerance = existing.uncertainty + uncertainty",
      "let maxMergeDistance = 0.10"],
     ["let tolerance = existing.uncertainty + uncertainty + SIMD3(repeating: 0.5)",
      "let maxMergeDistance = 100.0"],
     "CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM"),
    ("centring translation skipped in the expansion", SRC + "CIFImport.swift",
     "let mapped = wrapFractional(op.apply(site.position))",
     "let mapped = wrapFractional(op.rotation.apply(site.position))",
     "CIFImportTests.testBetaDoublePrimeImportsAs22AtomsC2OverM"),
    ("1/2+x read as x", SRC + "SymmetryOperation.swift",
     "constant += sign * magnitude",
     "constant += 0 * magnitude",
     "SymmetryOperationTests.testScrewAxisParsesRotationAndTranslation"),
    ("Laue table: 8-with-4-fold branch inverted", SRC + "PointSymmetry.swift",
     "case 8: return has4 ? .fourOverM : .mmm",
     "case 8: return has4 ? .mmm : .fourOverM",
     "PointSymmetryTests.testLaueClassesOfTheSmallerGroups"),
    ("detection without the target-distinct guard", SRC + "SymmetryDetection.swift",
     "!claimed.contains(index)\n                    && sites[index].z == site.z",
     "sites[index].z == site.z",
     "SymmetryDetectionTests.testSplitSiteGoldIsNotCubic"),
    ("detected translations not snapped to k/24", SRC + "SymmetryDetection.swift",
     "return abs(x - grid) <= tolerance ? grid : x",
     "return x",
     "PinsTests.testEveryStructureAgreesWithPymatgen"),
    ("the file's loop taken as the operations, not detected", SRC + "CIFImport.swift",
     "cell: cell, sites: parsed.sites, fileOperations: parsed.operations,",
     "cell: cell, sites: parsed.sites, fileOperations: nil,",
     "CIFImportTests.testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m"),
    ("proper-rotation orders read over every element", SRC + "PointSymmetry.swift",
     "let orders = Set(properRotations.compactMap(\\.order))",
     "let orders = Set(rotations.compactMap(\\.order))",
     "CIFImportTests.testTrigonalInHexagonalAxesIsThreeBarM"),
]


def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, **kw)


def main():
    if os.path.exists(SCRATCH):
        shutil.rmtree(SCRATCH)
    os.makedirs(SCRATCH)
    head = run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"]).stdout.strip()
    archive = subprocess.run(["git", "-C", ROOT, "archive", "HEAD"], stdout=subprocess.PIPE, check=True)
    subprocess.run(["tar", "-x", "-C", SCRATCH], input=archive.stdout, check=True)
    pristine = {}
    for _, file, old, _, _ in MUTANTS:
        path = os.path.join(SCRATCH, file)
        text = open(path).read()
        pristine[file] = text
        for piece in (old if isinstance(old, list) else [old]):
            count = text.count(piece)
            if count != 1:
                print(f"ABORT: the text for the mutant in {file} occurs {count} times, not once: {piece[:60]}")
                return 2
    print(f"mutants: {len(MUTANTS)} on commit {head}, scratch {SCRATCH}")
    # One warm build first, so each mutant run is an incremental rebuild.
    base = run(["swift", "test", "--filter", "DiffractoCoreTests"], cwd=SCRATCH)
    if base.returncode != 0:
        print("ABORT: the unmutated copy is not green")
        print(base.stdout[-3000:])
        return 2
    print("baseline: green")
    red = 0
    for name, file, old, new, test in MUTANTS:
        path = os.path.join(SCRATCH, file)
        mutated = pristine[file]
        for o, n in zip(old if isinstance(old, list) else [old], new if isinstance(new, list) else [new]):
            mutated = mutated.replace(o, n)
        open(path, "w").write(mutated)
        t0 = time.time()
        result = run(["swift", "test", "--filter", "DiffractoCoreTests"], cwd=SCRATCH)
        open(path, "w").write(pristine[file])
        failed = result.returncode != 0
        red += failed
        aimed = test.split(".")[-1]
        caught_by_aim = f"{aimed}]' failed" in result.stdout or f"error: -[DiffractoCoreTests.{test}]" in result.stdout
        verdict = "RED" if failed else "GREEN (test missing)"
        print(f"{verdict:20s} {name} [{file.split('/')[-1]}; aimed at {test}; "
              f"{'caught by it' if caught_by_aim else 'caught elsewhere' if failed else 'not caught'}; "
              f"{time.time() - t0:.0f} s]")
    print(f"mutants: {red} of {len(MUTANTS)} red on commit {head}")
    return 0 if red == len(MUTANTS) else 1


if __name__ == "__main__":
    sys.exit(main())
