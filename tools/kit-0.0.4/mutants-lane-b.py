#!/usr/bin/env python3
"""DiffractoKit 0.0.4 break-once run for G2 (registration section 3, lane B): every claim of RingCalibrationQTests is broken once.

Works in a COPY of the repository (never in the working tree): copies Package.swift, Sources and Tests to <copy>, applies ONE
source mutation at a time (the replaced text must occur exactly once), runs the filtered tests with
`swift test -c release --scratch-path <build> --filter <regex>`, records red or green, restores the file. A mutant that leaves
the filtered tests green is a hole in the tests. No sleeps, no polling. (Modelled on tools/slice-3b/mutants.py.)

  python3 tools/kit-0.0.4/mutants-lane-b.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
FILE = "Sources/DiffractoGeometry/RingCalibration+Q.swift"
FILTER = "RingCalibrationQTests"

# id, old text (exactly once), new text, what it breaks
MUTANTS = [
    ("factor-2-dropped",
     "let q = 2 * sin(atan2(radius, cameraLengthPixels) / 2) / wavelengthAngstroms",
     "let q = sin(atan2(radius, cameraLengthPixels) / 2) / wavelengthAngstroms",
     "q = sin(...)/lambda: the 2 of q = 2 sin(theta)/lambda dropped"),
    ("half-angle-dropped",
     "let q = 2 * sin(atan2(radius, cameraLengthPixels) / 2) / wavelengthAngstroms",
     "let q = 2 * sin(atan2(radius, cameraLengthPixels)) / wavelengthAngstroms",
     "q = 2 sin(atan(R/L))/lambda: the half angle theta = 2theta/2 dropped"),
    ("bragg-guard-removed",
     "        guard sinTheta < maximumLambdaOver2d else { return .nan }\n",
     "",
     "the inverse's q lambda/2 < 1/sqrt(2) guard removed"),
    ("overflow-guard-removed",
     "        return radius.isFinite && radius != 0 ? radius : .nan\n",
     "        return radius\n",
     "the inverse's NaN-on-overflow (and underflow) guard removed"),
    ("q-overflow-guard-removed",
     "        return q.isFinite && q > 0 ? q : .nan\n",
     "        return q\n",
     "the forward's NaN-on-overflow (and underflow) guard removed"),
    ("inverse-asin-without-factor-2",
     "let radius = cameraLengthPixels * tan(2 * asin(sinTheta))",
     "let radius = cameraLengthPixels * tan(asin(sinTheta))",
     "R = L tan(asin(q lambda/2)): the 2 of 2 asin dropped"),
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
    test = ["swift", "test", "-c", "release", "--scratch-path", str(build)]
    path = copy / FILE
    # the unmutated copy must be green, or a red mutant proves nothing
    base = run(test + ["--filter", FILTER], copy)
    print("baseline (unmutated copy): %s" % ("GREEN" if base.returncode == 0 else "RED (stop)"), flush=True)
    if base.returncode != 0:
        sys.exit(1)
    results = []
    for ident, old, new, what in MUTANTS:
        if wanted and ident not in wanted:
            continue
        original = path.read_text()
        if original.count(old) != 1:
            print("%s: the text to mutate occurs %d times, not once" % (ident, original.count(old)))
            results.append((ident, "NOT APPLIED"))
            continue
        path.write_text(original.replace(old, new))
        try:
            out = run(test + ["--filter", FILTER], copy)
        finally:
            path.write_text(original)
        text = out.stdout + out.stderr
        if out.returncode != 0 and "Test run" not in text:
            verdict = "BUILD ERROR"
        else:
            failed = sorted(set(re.findall(r"✘ Test (\w+)\(\)", text)))
            verdict = "RED (%d tests: %s)" % (len(failed), ", ".join(failed)) if out.returncode != 0 else "GREEN (hole)"
        print("%s [%s]: %s" % (ident, what, verdict), flush=True)
        results.append((ident, verdict))
    print("\nbreak-once: %d of %d mutants turned a test red" % (sum(1 for _, v in results if v.startswith("RED")), len(results)))


if __name__ == "__main__":
    main()
