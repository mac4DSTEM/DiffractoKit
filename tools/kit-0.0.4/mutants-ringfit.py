#!/usr/bin/env python3
"""DiffractoKit 0.0.4 break-once run for the ring-fit amendment, G1 and G7 (docs/registration-0.0.4.md §1.8, §2, §4): every
new test is broken once.

Works in a COPY of the repository (never in the working tree): copies Package.swift, Sources and Tests to <copy>, applies ONE
source mutation at a time (the replaced text must occur exactly once), runs the tests that must turn red with
`swift test -c release --scratch-path <build> --filter <regex>`, records red or green, restores the file. A mutant that leaves
the filtered tests green is a hole in the tests. Release (the ring fits are slow at -Onone); no sleeps, no polling.

  python3 tools/kit-0.0.4/mutants-ringfit.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
GEO = "Sources/DiffractoGeometry/"
ELL = GEO + "EllipseCalibration.swift"
EVI = GEO + "RingEvidence.swift"
CAL = GEO + "RingCalibration.swift"
FIT = GEO + "RingCalibration+Fit.swift"
AMEND = "RingFitAmendmentTests"
SPLIT = "RingFitSplitAndCancellationTests"

# id, file, old text (exactly once), new text, filter of the tests that must go red, what it breaks
MUTANTS = [
    ("evidence-coverage-dropped", CAL,
     "            || evidenceBins < RingEvidence.minimumOccupiedBins {", "            || false {",
     AMEND + "/aNoisyRingWithAWideGap", "the coverage refusal from the ring evidence removed (F1)"),
    ("occupancy-fraction", EVI,
     "public static let minimumAmplitudeFraction = 0.5", "public static let minimumAmplitudeFraction = 0.05",
     AMEND + "/aRingWithAWeakSector", "a sector holds the ring at 5 % of A instead of 50 %"),
    ("significance-raised", EVI,
     "public static let minimumSignificance = 3.0", "public static let minimumSignificance = 30.0",
     AMEND + "/theSameNoisyRingWithoutAGap", "a sector needs 30 standard errors instead of 3"),
    ("inner-drop-dropped", EVI,
     "return innerDrop >= Self.minimumDrop && outerDrop >= Self.minimumDrop",
     "return outerDrop >= Self.minimumDrop || innerDrop.isNaN",
     AMEND + "/aDiskEdge", "the peak test's inner drop not required (F2, disk edge)"),
    ("annulus-and-peak-dropped", CAL,
     "           !(fit.radius > inner && fit.radius < outer) {", "           false {",
     AMEND + "/theAnnulusRefusalReadsTheReportsRadii", "the annulus refusal removed (F2, blob; the peak test also refuses the blob)"),
    ("ellipse-shift-unbounded", EVI,
     "return !exceeds(centreShift, centreShiftError) && !exceeds(ellipseShift, ellipseShiftError)",
     "return !exceeds(centreShift, centreShiftError)",
     AMEND + "/aTexturedCircle", "the geometry check's ellipse-shift bound removed (F2, texture)"),
    ("impulse-removed", ELL,
     "        guard threshold > 0 else { return (found, threshold) }",
     "        guard threshold > 1e300 else { return (found, threshold) }",
     AMEND + "/anIsolatedImpulse", "the impulse exclusion removed (F3)"),
    ("impulse-factor", ELL,
     "return range.isFinite && range > 0 ? 2 * range : 0", "return range.isFinite && range > 0 ? 2000 * range : 0",
     AMEND + "/anIsolatedImpulse", "an impulse must exceed 2000 robust ranges instead of 2 (F3)"),
    ("normalisation-removed", ELL,
     "dynamicRange < 0.5 ? Double(sign: .plus, exponent: -dynamicRange.exponent - 1, significand: 1) : 1", "1",
     AMEND + "/theFitDoesNotDependOnScale", "no power-of-two normalisation (F6)"),
    ("signal-absolute", ELL,
     "&& dynamicRange > abs(maximum) * 1e-6", "&& dynamicRange > max(1, abs(maximum)) * 1e-6",
     AMEND + "/theFitDoesNotDependOnScale", "the signal test absolute again (F6)"),
    ("one-sided-floor-refused", ELL,
     "if atFloor == 2 || (atFloor == 1 && !requireRingEvidence) {", "if atFloor >= 1 {",
     AMEND + "/aRingWithOneSharpEdge", "one ring sigma at the floor still refused (F7)"),
    ("lowered-floor-dropped", ELL,
     "                ringSigmaFloor = 0.01\n", "                ringSigmaFloor = 0.25\n",
     AMEND + "/aRingWithOneSharpEdge", "the one-sided refit keeps the 0.25 px floor (F7's radius bias)"),
    ("noHighTension-not-appended", FIT,
     "if highTensionVolts == nil, !invalid { refused.append(.noHighTension) }", "_ = invalid",
     SPLIT + "/theSplitEqualsTheFit", "calibrating drops .noHighTension (G1)"),
    ("fitRing-without-evidence", FIT,
     "requireRingEvidence: true, progress: progress, isCancelled: isCancelled)",
     "requireRingEvidence: false, progress: progress, isCancelled: isCancelled)",
     SPLIT + "/theSplitEqualsTheFit|" + AMEND + "/aNoisyRingWithAWideGap", "fitRing without the ring evidence (G1, §1.4)"),
    ("cancellation-ignored", ELL,
     "if Task.isCancelled || isCancelled?() == true { throw CancellationError() }", "_ = isCancelled",
     SPLIT + "/progressIsMonotone", "the cancellation check removed (G7)"),
    ("final-progress-dropped", ELL,
     "                                    innerRadius: innerRadius, outerRadius: outerRadius)\n        span(1)\n",
     "                                    innerRadius: innerRadius, outerRadius: outerRadius)\n",
     SPLIT + "/progressIsMonotone", "the final progress report dropped (G7)"),
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
    test = ["swift", "test", "-c", "release", "--scratch-path", str(build)]
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
            failed = sorted(set(re.findall(r"✘ Test (\w+)\(", text)))
            verdict = "RED (%d tests: %s)" % (len(failed), ", ".join(failed)) if out.returncode != 0 else "GREEN (hole)"
        print("%s [%s]: %s" % (ident, what, verdict), flush=True)
        results.append((ident, verdict))
    print("\nbreak-once: %d of %d mutants turned a test red" % (sum(1 for _, v in results if v.startswith("RED")), len(results)))


if __name__ == "__main__":
    main()
