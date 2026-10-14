#!/usr/bin/env python3
"""Break-once proof for slice 3's Gate B fixes (registration-slice-3 section 7; the lanes' own mutants ran in their copies, logs in the session scratchpad): each planned mutant is applied to a COPY of the
committed tree (git archive HEAD into a scratch directory; the repository is never edited), DiffractoCoreTests is run,
and the run must go RED. A mutant that stays green means a test is missing. Prints one line per mutant and a summary;
exits non-zero unless every mutant is red.

    python3 -I tools/slice-3/mutants.py [<scratch-dir>]        (default: <repo>/.build/mutants)

Each mutant is (name, file, exact old text, new text, the test it is aimed at); the old text must occur exactly once.
"""
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRATCH = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, ".build", "mutants")
SRC = "Sources/DiffractoGeometry/"

MUTANTS = [
    ("F1: no overflow guard on the wavelength", SRC + "ElectronWavelength.swift",
     "guard product.isFinite else { return .nan }",
     "_ = 0",
     "GeometryRefuterTests.testWavelengthOverflowIsNaN"),
    ("F2: camera length may be infinite", SRC + "RingCalibration.swift",
     "return length.isFinite ? length : .nan",
     "return length",
     "GeometryRefuterTests.testOverflowingLengthsAreNaN"),
    ("F2: dq may be infinite", SRC + "RingCalibration.swift",
     "return q.isFinite ? q : .nan",
     "return q",
     "GeometryRefuterTests.testOverflowingLengthsAreNaN"),
    ("F4: project accepts a non-unit n", SRC + "ZoneAxis.swift",
     "guard isUnitDirection(n) else { return [] }",
     "_ = 0",
     "GeometryRefuterTests.testProjectRefusesBadDirectionsAndWavelengths"),
    ("F5: project accepts a negative wavelength", SRC + "ZoneAxis.swift",
     "if let lambda = wavelengthAngstroms { guard lambda.isFinite, lambda >= 0 else { return [] } }",
     "_ = 0",
     "GeometryRefuterTests.testProjectRefusesBadDirectionsAndWavelengths"),
    ("F7: the list overload with a free k_max of 4", SRC + "RingIndexing.swift",
     "try index(q: q, sigma: sigma, in: list.families, kMax: list.kMax, coverage: k)",
     "try index(q: q, sigma: sigma, in: list.families, kMax: 4, coverage: k)",
     "GeometryRefuterTests.testListOverloadUsesTheListsKMax"),
]


def touch_module():
    """Every source of DiffractoCore gets a new mtime, so the debug build recompiles the whole module. Without it an
    incremental build recompiles only the mutated file, and a constant another file folded in (a default argument
    reading `defaultMetricLengthTolerance`) keeps its old value: the first run of 2026-10-15 reported such mutants
    green although the code they change was never rebuilt."""
    now = time.time()
    folder = os.path.join(SCRATCH, SRC)
    for name in os.listdir(folder):
        os.utime(os.path.join(folder, name), (now, now))


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
    base = run(["swift", "test", "--filter", "DiffractoGeometryTests"], cwd=SCRATCH)
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
        touch_module()
        t0 = time.time()
        try:
            result = run(["swift", "test", "--filter", "DiffractoGeometryTests"], cwd=SCRATCH, timeout=900)
        except subprocess.TimeoutExpired as e:   # a mutant that hangs the tests counts as caught, as a timeout
            result = subprocess.CompletedProcess(e.cmd, 124, stdout="(timeout)")
        open(path, "w").write(pristine[file])
        touch_module()
        failed = result.returncode != 0
        red += failed
        aimed = test.split(".")[-1]
        caught_by_aim = (f"{aimed}]' failed" in result.stdout or f"error: -[DiffractoCoreTests.{test}]" in result.stdout
                         or f"Test {aimed}() failed" in result.stdout or f"Test {aimed}(" in result.stdout and "failed" in result.stdout
                         and any(line.startswith("✘") and f"{aimed}(" in line for line in result.stdout.splitlines()))
        verdict = "RED" if failed else "GREEN (test missing)"
        print(f"{verdict:20s} {name} [{file.split('/')[-1]}; aimed at {test}; "
              f"{'caught by it' if caught_by_aim else 'caught elsewhere' if failed else 'not caught'}; "
              f"{time.time() - t0:.0f} s]")
    print(f"mutants: {red} of {len(MUTANTS)} red on commit {head}")
    return 0 if red == len(MUTANTS) else 1


if __name__ == "__main__":
    sys.exit(main())
