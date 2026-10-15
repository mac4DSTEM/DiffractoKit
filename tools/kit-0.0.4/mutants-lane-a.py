#!/usr/bin/env python3
"""DiffractoKit 0.0.4 break-once run, lane A (tetragonal along a or b, section 6; error texts, section 5).

Works in a COPY of the repository (never in the working tree): copies Package.swift, Sources and Tests to <copy>, runs the
filtered tests once unmutated (they must be green), then applies ONE source mutation at a time (the replaced text must occur
exactly once), builds and runs the tests that must turn red with `swift test -c release --scratch-path <build> --filter <regex>`,
records red or green, restores the file. A mutant that leaves the filtered tests green is a hole in the tests; a compile error
is reported as BUILD ERROR and is not red. No sleeps, no polling. (`swift build --build-tests -c release` is never used.)

  python3 tools/kit-0.0.4/mutants-lane-a.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
UC = "Sources/DiffractoCore/UnitCell.swift"
PS = "Sources/DiffractoCore/PointSymmetry.swift"
EXT = "public var errorDescription: String? { description }"
CIF = "Sources/DiffractoCore/CIFImport.swift"
RL = "Sources/DiffractoCore/ReflectionList.swift"
KE = "Sources/DiffractoScattering/StructureFactors.swift"

TETRA = "TetragonalSettingTests|UnitCellTests|PointSymmetryTests"
FAMILY_B = "        if ac && !ab && ra && rb && rc { return .tetragonal(uniqueAxis: .b) }\n"
FAMILY_A = "        if bc && !ab && ra && rb && rc { return .tetragonal(uniqueAxis: .a) }\n"

# id, file, old text (exactly once), new text, filter of the tests that must go red, what it breaks
MUTANTS = [
    ("family-c-only", UC, FAMILY_B + FAMILY_A, "",
     TETRA, "the tetragonal fix reverted: only a = b is tetragonal (4-fold along c)"),
    ("family-b-dropped", UC, FAMILY_B, "",
     TETRA, "a = c != b (4-fold along b) read orthorhombic again"),
    ("family-a-dropped", UC, FAMILY_A, "",
     TETRA, "b = c != a (4-fold along a) read orthorhombic again"),
    ("holohedry-a-about-c", PS, "case .a: generators = [c4x, c2y, .inversion]", "case .a: generators = [c4z, c2x, .inversion]",
     TETRA, "the holohedry of a tetragonal-along-a metric built about c"),
    ("holohedry-b-about-c", PS, "case .b: generators = [c4y, c2z, .inversion]", "case .b: generators = [c4z, c2x, .inversion]",
     TETRA, "the holohedry of a tetragonal-along-b metric built about c"),
    ("holohedry-a-about-b", PS, "case .a: generators = [c4x, c2y, .inversion]", "case .a: generators = [c4y, c2z, .inversion]",
     TETRA, "the holohedry along a built about b (wrong, but still a 4/mmm: the isometry and the axis tests catch it)"),
    ("cif-error-nil", CIF, EXT, "public var errorDescription: String? { nil }",
     "ErrorTextTests", "CIFImportError.errorDescription returns nil"),
    ("cif-error-fixed", CIF, EXT, 'public var errorDescription: String? { "CIF import failed." }',
     "ErrorTextTests", "CIFImportError.errorDescription returns a fixed string"),
    ("reflection-error-nil", RL, EXT, "public var errorDescription: String? { nil }",
     "ErrorTextTests", "ReflectionListError.errorDescription returns nil"),
    ("reflection-error-fixed", RL, EXT, 'public var errorDescription: String? { "Reflection list failed." }',
     "ErrorTextTests", "ReflectionListError.errorDescription returns a fixed string"),
    ("kinematic-error-nil", KE, EXT, "public var errorDescription: String? { nil }",
     "KinematicErrorTextTests", "KinematicError.errorDescription returns nil"),
    ("kinematic-error-fixed", KE, EXT, 'public var errorDescription: String? { "Kinematic calculation failed." }',
     "KinematicErrorTextTests", "KinematicError.errorDescription returns a fixed string"),
]
BASELINE_FILTER = TETRA + "|ErrorTextTests|KinematicErrorTextTests"


def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def verdict_of(out):
    text = out.stdout + out.stderr
    failed = sorted(set(re.findall(r"✘ Test (\w+)\(\)", text)))
    ran = "Test run with" in text
    if out.returncode == 0 and ran:
        return "GREEN"
    if failed:
        return "RED (%d tests: %s)" % (len(failed), ", ".join(failed))
    return "BUILD ERROR"


def main():
    copy, build = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    wanted = set(sys.argv[3:])
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    shutil.copy(ROOT / "Package.swift", copy / "Package.swift")
    for d in ("Sources", "Tests"):
        shutil.copytree(ROOT / d, copy / d)
    test = ["swift", "test", "-c", "release", "-j", "8", "--scratch-path", str(build)]
    base = verdict_of(run(test + ["--filter", BASELINE_FILTER], copy))
    print("baseline (unmutated copy): %s" % base, flush=True)
    if base != "GREEN":
        sys.exit("the unmutated copy is not green: stop")
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
        verdict = verdict_of(out)
        if verdict == "GREEN":
            verdict = "GREEN (hole)"
        print("%s [%s]: %s" % (ident, what, verdict), flush=True)
        results.append((ident, verdict))
    print("\nbreak-once: %d of %d mutants turned a test red" % (sum(1 for _, v in results if v.startswith("RED")), len(results)))


if __name__ == "__main__":
    main()
