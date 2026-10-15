#!/usr/bin/env python3
"""DiffractoKit 0.2.0 break-once run for the policy parameters (docs/registration-0.2.0.md §2): one mutant per policy group,
each changing the group's DEFAULT; the full suite must go red.

Works in a COPY of the repository (never in the working tree): copies Package.swift, Sources, Tests and Examples to <copy>,
applies ONE source mutation at a time (the replaced text must occur exactly once), runs `swift test -c release` there with
its own scratch path, records red or green, restores the file. A mutant that leaves the suite green is a hole in the tests.
A verdict says which tests turned red and how many of them are the policy tests' own default-value checks (the `theDefaults...`
tests of Tests/*/PolicyTests.swift): a mutant caught ONLY by those is caught by a registered default, not by a behaviour test,
and is reported as such. No sleeps, no polling.

  python3 tools/kit-0.2.0/mutants-policy.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
CORE = "Sources/DiffractoCore/"
SCAT = "Sources/DiffractoScattering/"
GEO = "Sources/DiffractoGeometry/"
EPOL = GEO + "EllipseCalibration+Policy.swift"
RCAL = GEO + "RingCalibration.swift"
REVI = GEO + "RingEvidence.swift"
RIDX = GEO + "RingIndexing.swift"
RLST = CORE + "ReflectionList.swift"
CSTR = CORE + "CrystalStructure.swift"
SDET = CORE + "SymmetryDetection.swift"
SOPR = CORE + "SymmetryOperation.swift"
CIFI = CORE + "CIFImport.swift"
SFAC = SCAT + "StructureFactors.swift"

# id, group, file, old text (exactly once), new text, what it breaks
MUTANTS = [
    ("ellipse-coverage", "ring fit: coverage", EPOL,
     "degeneracyBoundBins ?? angularBinCount * 5 / 6", "degeneracyBoundBins ?? angularBinCount * 2 / 3",
     "the full-coverage bound 30 of 36 becomes 24"),
    ("ellipse-one-ring", "ring fit: one-ring ratio", EPOL,
     "oneRingRatio: Double = 1.10,", "oneRingRatio: Double = 1.5,", "the one-ring ratio 1.10 becomes 1.5"),
    ("ellipse-profile-residual", "ring fit: residual caps", EPOL,
     "maximumProfileResidual: Double = 0.22,", "maximumProfileResidual: Double = 0.30,", "the profile residual cap 0.22 becomes 0.30"),
    ("ellipse-floor-tolerance", "ring fit: ring-width floor", EPOL,
     "ringWidthFloorTolerance: Double = 0.05,", "ringWidthFloorTolerance: Double = 0.0,", "the 0.30 px guard becomes the 0.25 px floor"),
    ("ellipse-robust-cut", "ring fit: robust refit", EPOL,
     "robustCutMultiplier: Double = 4 * 1.4826,", "robustCutMultiplier: Double = 1000,", "the robust refit never finds an outlier"),
    ("ellipse-impulse", "ring fit: impulse exclusion", EPOL,
     "impulseRangeFactor: Double = 2,", "impulseRangeFactor: Double = 2000,", "an impulse must exceed 2000 robust ranges"),
    ("ellipse-axis-ratio", "ring fit: axis ratio", EPOL,
     "maximumAxisRatio: Double = 5,", "maximumAxisRatio: Double = 50,", "the axis ratio 5 becomes 50"),
    ("ellipse-signal", "ring fit: relative signal", EPOL,
     "minimumRelativeSignal: Double = 1e-6,", "minimumRelativeSignal: Double = 1e-1,", "the relative signal test 1e-6 becomes 0.1"),
    ("calibration-residual", "calibration: refusal bounds", RCAL,
     "public static let maximumResidual = 0.22", "public static let maximumResidual = 0.5", "the calibration residual cap 0.22 becomes 0.5"),
    ("calibration-ellipticity", "calibration: refusal bounds", RCAL,
     "public static let maximumEllipticity = 0.093", "public static let maximumEllipticity = 0.5", "the ellipticity bound 0.093 becomes 0.5"),
    ("calibration-coverage", "calibration: refusal bounds", RCAL,
     "public static let minimumOccupiedAngularBins = 30", "public static let minimumOccupiedAngularBins = 20",
     "the calibration's occupied-bin bound 30 becomes 20"),
    ("evidence-drop", "ring evidence", REVI,
     "public static let minimumDrop = 0.25", "public static let minimumDrop = 0.0", "the peak test's drop 0.25 becomes 0"),
    ("evidence-fraction", "ring evidence", REVI,
     "public static let minimumAmplitudeFraction = 0.5", "public static let minimumAmplitudeFraction = 0.05",
     "a sector holds the ring at 5 % of the amplitude"),
    ("indexing-margin", "ring indexing", RIDX,
     "public static let minimumMargin: Double = log(19.0)", "public static let minimumMargin: Double = log(2.0)", "the margin ln 19 becomes ln 2"),
    ("list-kmax", "reflection list: refusal bounds", RLST,
     "public static let maximumKMax: Double = 4", "public static let maximumKMax: Double = 5", "the k_max cap 4 becomes 5"),
    ("list-reflections", "reflection list: refusal bounds", RLST,
     "public static let maximumReflections = 10_000_000", "public static let maximumReflections = 1_000",
     "the reflection cap 10 million becomes 1000"),
    ("list-family-tolerance", "reflection list: ring rule", RLST,
     "public static let familyTolerance: Double = 1e-9", "public static let familyTolerance: Double = 1e-4",
     "the ring-grouping tolerance 1e-9 becomes 1e-4"),
    ("list-absence-tolerance", "reflection list: absence rule", RLST,
     "absenceTolerance: Double = 1e-9)", "absenceTolerance: Double = 1e-2)", "the absence tolerance 1e-9 becomes 1e-2"),
    ("validation-contact", "structure validation", CSTR,
     "public static let closeContactLimit = 0.5", "public static let closeContactLimit = 0.05", "the close-contact limit 0.5 becomes 0.05"),
    ("validation-split", "structure validation", CSTR,
     "splitSiteOccupancyTolerance: Double = 0.02)", "splitSiteOccupancyTolerance: Double = 0.5)", "the split-site sum 1.02 becomes 1.5"),
    ("detect-occupancy", "symmetry tolerances", SDET,
     "public static let occupancyTolerance = 0.02", "public static let occupancyTolerance = 0.5", "the occupancy tolerance 0.02 becomes 0.5"),
    ("detect-grid", "symmetry tolerances", SDET,
     "public static let translationGrid = 1.0 / 24", "public static let translationGrid = 1.0 / 5", "the translation grid 1/24 becomes 1/5"),
    ("detect-precision", "symmetry tolerances", SDET,
     "gridSnapFactor: Double = 3, precisionLimitFactor: Double = 10", "gridSnapFactor: Double = 3, precisionLimitFactor: Double = 0",
     "the precision-limited factor 10 becomes 0 (nothing reported unresolved)"),
    ("operation-magnitude", "CIF", SOPR,
     "defaultMaximumMagnitude: Double = 1000", "defaultMaximumMagnitude: Double = 1e12", "the coefficient bound 1000 becomes 1e12"),
    ("cif-merge-distance", "CIF", CIFI,
     "maximumMergeDistance: Double = 0.10,", "maximumMergeDistance: Double = 0.001,", "the merge distance 0.10 A becomes 0.001"),
    ("cif-half-step-cap", "CIF", CIFI,
     "maximumCoordinateHalfStep: Double = 1e-2,", "maximumCoordinateHalfStep: Double = 1e-4,", "the half-step cap 1e-2 becomes 1e-4"),
    ("cif-right-angle", "CIF", CIFI,
     "rightAngleToleranceDeg: Double = 0.05,\n                    maximumOperationMagnitude", "rightAngleToleranceDeg: Double = 0.0,\n                    maximumOperationMagnitude",
     "the right-angle tolerance 0.05 becomes 0"),
    ("scattering-cut", "|F| cut", SFAC,
     "public static let defaultCut = 1e-4", "public static let defaultCut = 1e-2", "the |F| cut 1e-4 becomes 1e-2"),
]


def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def main():
    copy, build = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    wanted = set(sys.argv[3:])
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    shutil.copy(ROOT / "Package.swift", copy / "Package.swift")
    for d in ("Sources", "Tests", "Examples"):
        shutil.copytree(ROOT / d, copy / d)
    own = set()
    for f in (ROOT / "Tests").glob("*/PolicyTests.swift"):
        own |= set(re.findall(r"func (\w+)\(", f.read_text()))
    test = ["swift", "test", "-c", "release", "--scratch-path", str(build), "-j", "2"]
    results = []
    for ident, group, rel, old, new, what in MUTANTS:
        if wanted and ident not in wanted:
            continue
        path = copy / rel
        original = path.read_text()
        if original.count(old) != 1:
            print("%s: the text to mutate occurs %d times, not once" % (ident, original.count(old)))
            results.append((ident, group, "NOT APPLIED"))
            continue
        path.write_text(original.replace(old, new))
        try:
            out = run(test, copy)
        finally:
            path.write_text(original)
        text = out.stdout + out.stderr
        failed = sorted(set(re.findall(r"✘ Test (\w+)\(", text)))
        if out.returncode == 0:
            verdict = "GREEN (hole)"
        elif not failed:
            verdict = "BUILD ERROR"
        else:
            policy_only = [n for n in failed if n in own]
            behaviour = [n for n in failed if n not in own]
            verdict = "RED (%d tests; %d behaviour tests: %s; %d policy-file tests: %s)" % (
                len(failed), len(behaviour), ", ".join(behaviour[:6]) + (" ..." if len(behaviour) > 6 else ""),
                len(policy_only), ", ".join(policy_only[:6]))
        print("%s [%s; %s]: %s" % (ident, group, what, verdict), flush=True)
        results.append((ident, group, verdict))
    print("\nbreak-once: %d of %d mutants turned a test red" % (sum(1 for _, _, v in results if v.startswith("RED")), len(results)))
    groups = sorted({g for _, g, _ in results})
    for g in groups:
        reds = [v for _, gg, v in results if gg == g and v.startswith("RED")]
        print("  group %-36s %d of %d red" % (g, len(reds), sum(1 for _, gg, _ in results if gg == g)))


if __name__ == "__main__":
    main()
