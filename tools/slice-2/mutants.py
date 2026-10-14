#!/usr/bin/env python3
"""Break-once proof for slice 2 (registration-slice-2 section 6, gate 3; slice 1's script, its mutant list replaced): each planned mutant is applied to a COPY of the
committed tree (git archive HEAD into a scratch directory; the repository is never edited), DiffractoCoreTests is run,
and the run must go RED. A mutant that stays green means a test is missing. Prints one line per mutant and a summary;
exits non-zero unless every mutant is red.

    python3 -I tools/slice-2/mutants.py [<scratch-dir>]        (default: <repo>/.build/mutants)

Each mutant is (name, file, exact old text, new text, the test it is aimed at); the old text must occur exactly once.
"""
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRATCH = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, ".build", "mutants")
SRC = "Sources/DiffractoScattering/"

MUTANTS = [
    ("table row Z-1 -> Z", SRC + "ElectronScatteringFactors.swift",
     "let p = parameters[z - 1]",
     "let p = parameters[min(z, maximumAtomicNumber - 1)]",
     "ScatteringPinsTests.testStructureMatchesNumPyPins"),
    ("f_e of |g| instead of |g|^2", SRC + "StructureFactors.swift",
     "let gSq = gLen * gLen",
     "let gSq = gLen",
     "ScatteringTruthTests.testFCCExtinctionsAndClosedForm"),
    ("no 1/V", SRC + "StructureFactors.swift",
     "re /= volume; im /= volume",
     "_ = volume",
     "ScatteringTruthTests.testFCCExtinctionsAndClosedForm"),
    ("cut >= instead of >", SRC + "StructureFactors.swift",
     "let kept = all.filter { $0.structureFactor.magnitude > cut }",
     "let kept = all.filter { $0.structureFactor.magnitude >= cut }",
     "ScatteringPinsTests.testCutEdgeKeepsNothing"),
    ("default cut 1e-3", SRC + "StructureFactors.swift",
     "public static let defaultCut = 1e-4",
     "public static let defaultCut = 1e-3",
     "ScatteringTruthTests.testWeakCutStraddlesTheCut"),
    ("occupancy ignored", SRC + "StructureFactors.swift",
     "re += fe * site.occupancy * cos(phase)\n            im += fe * site.occupancy * sin(phase)",
     "re += fe * cos(phase)\n            im += fe * sin(phase)",
     "ScatteringTruthTests.testOccupancyClosedForm"),
    ("phase sign flipped", SRC + "StructureFactors.swift",
     "im += fe * site.occupancy * sin(phase)",
     "im -= fe * site.occupancy * sin(phase)",
     "ScatteringPinsTests.testZincBlendeImaginaryPartSign"),
    ("last site skipped", SRC + "StructureFactors.swift",
     "for site in sites {\n            guard let fe = formFactor(site.z, gSq) else { return nil }",
     "for site in sites.dropLast() {\n            guard let fe = formFactor(site.z, gSq) else { return nil }",
     "ScatteringTruthTests.testDiamondExtinctionsAndClosedForm"),
    ("families on the uncut list", SRC + "StructureFactors.swift",
     "let rings = ReflectionList.families(of: kept.map(\\.reflection), laueGroup: list.laueGroup,",
     "let rings = ReflectionList.families(of: all.map(\\.reflection), laueGroup: list.laueGroup,",
     "ScatteringTruthTests.testRingCountsAfterTheCut"),
    ("unsupported elements not checked first", SRC + "StructureFactors.swift",
     "guard missing.isEmpty else { throw .unsupportedElements(missing) }\n        let issues",
     "let issues",
     "ScatteringTruthTests.testRefusals"),
    ("Gate B F2: the structure's validation not called", SRC + "StructureFactors.swift",
     "guard issues.isEmpty else { throw .invalidStructure(issues.map(\\.message)) }",
     "_ = issues",
     "ScatteringRefuterTests.testInvalidStructuresAreRefused"),
    ("Gate B F2: non-finite F accepted", SRC + "StructureFactors.swift",
     "guard f.real.isFinite, f.imaginary.isFinite else { throw .nonFiniteStructureFactor(r.hkl) }",
     "_ = 0",
     "ScatteringRefuterTests.testNonFiniteFormFactorsAreRefused"),
    ("Gate B F4: keptButOperationAbsent filters nothing", SRC + "StructureFactors.swift",
     "public var keptButOperationAbsent: [KinematicReflection] { kept.filter { !$0.reflection.isAllowed } }",
     "public var keptButOperationAbsent: [KinematicReflection] { [] }",
     "ScatteringRefuterTests.testOperationAbsentReflectionsAboveTheCutAreNamed"),
    ("Gate B F8: f_e for a negative g2", SRC + "ElectronScatteringFactors.swift",
     "guard isSupported(z: z), gSquared.isFinite, gSquared >= 0 else { return nil }",
     "guard isSupported(z: z) else { return nil }",
     "ScatteringRefuterTests.testScatteringFactorDomain"),
    ("ring intensity not summed (first member only)", SRC + "StructureFactors.swift",
     "KinematicFamily(family: family, intensity: family.members.reduce(0) { $0 + (byHKL[$1] ?? 0) })",
     "KinematicFamily(family: family, intensity: byHKL[family.members[0]] ?? 0)",
     "ScatteringPinsTests.testStructureMatchesNumPyPins"),
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
    base = run(["swift", "test", "--filter", "DiffractoScatteringTests"], cwd=SCRATCH)
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
            result = run(["swift", "test", "--filter", "DiffractoScatteringTests"], cwd=SCRATCH, timeout=900)
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
