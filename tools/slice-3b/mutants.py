#!/usr/bin/env python3
"""DiffractoKit slice 3b break-once run (registration section 5, gate 3): every claim of the 3b tests is broken once.

Works in a COPY of the repository (never in the working tree): copies Package.swift, Sources and Tests to <copy>, applies ONE
source mutation at a time (the replaced text must occur exactly once), builds and runs the tests that must turn red with
`swift test --scratch-path <build> --filter <regex>`, records red or green, restores the file. A mutant that leaves the filtered
tests green is a hole in the tests. Debug build (the instruction for these lanes); no sleeps, no polling.

  python3 tools/slice-3b/mutants.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
GEO = "Sources/DiffractoGeometry/"

# id, file, old text (exactly once), new text, filter of the tests that must go red (the Debug parity tests take 15 minutes and are
# left out of these filters; the parity sensitivity is shown by tools/slice-3b/parity-mutants.py), what it breaks
MUTANTS = [
    ("coverage-constant", GEO + "EllipseCalibration.swift",
     "public static let degeneracyBoundBins = angularBinCount * 5 / 6", "public static let degeneracyBoundBins = angularBinCount * 2 / 3",
     "BinConstants|AGappedRingIsRefused|RingFitCalibrationTests", "the 30-bin coverage bound lowered to 24"),
    ("floor-0.30", GEO + "EllipseCalibration.swift",
     "min(parameters[3], parameters[4]) > 0.25 + 0.05", "min(parameters[3], parameters[4]) > 0.25",
     "ARingFittedWithinPointThreePx", "ADR 069's ring-width floor guard 0.30 px -> 0.25 px"),
    ("second-start-dropped", GEO + "EllipseCalibration.swift",
     "abs(ringSigma - parameters[3]) > 1e-9", "abs(ringSigma - parameters[3]) > 1e9",
     "NoisyGappedRing|NoisyRingDoesNotDrift", "ADR 069's half-maximum start never added"),
    ("profile-residual-cap", GEO + "EllipseCalibration.swift",
     "maximumProfileResidual: Double = 0.22,", "maximumProfileResidual: Double = 0.30,",
     "TheProfileResidualCap", "fitBestAvailable's profile residual cap 0.22 -> 0.30"),
    ("annulus-outer-boundary", GEO + "EllipseCalibration.swift",
     "guard radius > innerRadius, radius <= outerRadius else { continue }", "guard radius > innerRadius, radius < outerRadius else { continue }",
     "TheAnnulusBoundaries", "fit1D's annulus outer boundary <= -> <"),
    ("transform-e", GEO + "EllipseTransform.swift",
     "let e = ellipse.b / ellipse.a", "let e = ellipse.a / ellipse.b",
     "EllipseTransformTests", "the transform's e = b/a inverted"),
    ("stop-rule", GEO + "EllipseCalibration.swift",
     "                if relativeStep < 1e-10 {\n                    converged = true\n                    break\n                }\n                damping *= 10\n            }\n        }\n        guard converged",
     "                damping *= 10\n            }\n        }\n        guard converged",
     "EllipseStopRuleTests", "the stop test on a refused step removed (ADR-less Gate D of 2026-10-12)"),
    ("jacobian-b", GEO + "EllipseCalibration.swift",
     "                    dx * dy * value,\n", "                    dx * dy,\n",
     "CleanSweepEqualsSciPy", "fit1D's B column of the Jacobian unweighted (wrong gradient: SciPy's optimum is not reached)"),
    ("fit-radius-from-a", GEO + "RingCalibration+Fit.swift",
     "calibrate(radius: report.radius,", "calibrate(radius: report.semiAxisA,",
     "RingFitCalibrationTests", "the calibration's radius taken as the semi-axis a instead of (a + b)/2"),
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
    for d in ("Sources", "Tests"):
        shutil.copytree(ROOT / d, copy / d)
    # the copy's tests read fixtures by #filePath, which is inside the copy; the copy is a full copy
    test = ["swift", "test", "--scratch-path", str(build)]
    results = []
    for ident, rel, old, new, flt, what in MUTANTS:
        if wanted and ident not in wanted:
            continue
        path = copy / rel
        original = path.read_text()
        if original.count(old) != 1:
            print("%s: the text to mutate occurs %d times, not once" % (ident, original.count(old)))
            results.append((ident, "NOT APPLIED"))
            continue
        path.write_text(original.replace(old, new))
        try:
            out = run(test + ["--filter", flt], copy)
        finally:
            path.write_text(original)
        text = out.stdout + out.stderr
        if "error:" in text and "Compiling" not in text and "recorded an issue" not in text and out.returncode != 0 and "Test run" not in text:
            verdict = "BUILD ERROR"
        else:
            failed = sorted(set(re.findall(r"✘ Test (\w+)\(\)", text)))
            verdict = "RED (%d tests: %s)" % (len(failed), ", ".join(failed)) if out.returncode != 0 else "GREEN (hole)"
        print("%s [%s]: %s" % (ident, what, verdict), flush=True)
        results.append((ident, verdict))
    print("\nbreak-once: %d of %d mutants turned a test red" % (sum(1 for _, v in results if v.startswith("RED")), len(results)))


if __name__ == "__main__":
    main()
