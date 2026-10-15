#!/usr/bin/env python3
"""DiffractoKit slice 1b break-once run (registration docs/registration-slice-1b.md section 4; slice 1's pattern).

Works in a COPY (never in the working tree) that holds only what slice 1b touches: Package.swift reduced to DiffractoCore and
DiffractoCoreTests, Sources/DiffractoCore, Tests/DiffractoCoreTests, and slice 2's pins file (read by path by one test). The
new tests (ReductionTests, ConventionalCellTests, ConventionalSymmetryTests, FastDetectionTests) must be green unmutated;
then ONE source mutation at a time (the replaced text must occur exactly once) is applied, the new tests are run with
`swift test -c release -j 2 --scratch-path <build> --filter <regex>`, red or green is recorded and the file restored. A mutant
that leaves them green is a hole (or, when argued in the report, an equivalent mutant); a compile error is BUILD ERROR and is
not red. No sleeps, no polling. (`swift build --build-tests -c release` is never used.)

  python3 tools/slice-1b/mutants.py <copy-dir> <build-scratch-dir> [mutant-id ...]
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RED = "Sources/DiffractoCore/LatticeReduction.swift"
CON = "Sources/DiffractoCore/ConventionalCell.swift"
SET = "Sources/DiffractoCore/ConventionalSetting.swift"
DET = "Sources/DiffractoCore/SymmetryDetection.swift"
STR = "Sources/DiffractoCore/CrystalStructure.swift"

NEW = "ReductionTests|ConventionalCellTests|ConventionalSymmetryTests|FastDetectionTests"

# id, file, old text (exactly once), new text, what it breaks
MUTANTS = [
    # --- the Niggli reduction
    ("niggli-length-order-flipped", RED, "if a > b + eps || (abs(a - b) <= eps && abs(xi) > abs(eta) + eps) { cols.swapAt(0, 1); continue }",
     "if a < b - eps || (abs(a - b) <= eps && abs(xi) > abs(eta) + eps) { cols.swapAt(0, 1); continue }",
     "A1: the sign of the length condition A <= B flipped"),
    ("niggli-epsilon-x1e6", RED, "public static let niggliRelativeEpsilon = 1e-5", "public static let niggliRelativeEpsilon = 10.0",
     "epsilon 1e-5 -> 1e1 (x1e6): every comparison a tie"),
    ("niggli-type-II-sign", RED, "for i in 0..<3 { e[i] = s[i] == 0 ? 1 : -s[i] }", "for i in 0..<3 { e[i] = s[i] == 0 ? 1 : s[i] }",
     "A4: the all-non-positive form made from the wrong signs"),
    ("niggli-A8-dropped", RED, "if sum < -eps || (abs(sum) <= eps && 2 * (a + eta) + zeta > eps) {", "if false {",
     "A8: c -> a + b + c never taken"),
    ("niggli-A7-dropped", RED, "if abs(zeta) > a + eps ||", "if false && abs(zeta) > a + eps ||",
     "A7: the shear of b by a never taken"),
    ("niggli-shear-sign", RED, "cols[2] = cols[2] &- (xi >= 0 ? 1 : -1) &* cols[1]; continue", "cols[2] = cols[2] &+ (xi >= 0 ? 1 : -1) &* cols[1]; continue",
     "A5: the shear of c by b with the wrong sign"),
    ("niggli-handedness", RED, "if p.determinant < 0 { p = p.negated }", "",
     "the determinant -1 basis is not made right-handed"),
    # --- the conventional cell
    ("char-tI-tP-swapped", CON, "case 16: bravais = d == 1 ? .tP : .tI", "case 16: bravais = d == 1 ? .tI : .tP",
     "the tetragonal centring condition swapped"),
    ("char-cI-cF-swapped", CON, "case 48: bravais = d == 1 ? .cP : (d == 2 ? .cI : .cF)", "case 48: bravais = d == 1 ? .cP : (d == 2 ? .cF : .cI)",
     "the cubic centring condition swapped"),
    ("char-oI-oC-condition", CON, "classes.contains(SIMD3(1, 1, 1))", "classes.contains(SIMD3(1, 1, 0))",
     "the body-centring test of an orthorhombic cell reads a face"),
    ("char-hR-reverse", CON, "if !pointClasses(pc).classes.contains(SIMD3(2, 1, 1)) {", "if pointClasses(pc).classes.contains(SIMD3(2, 1, 1)) {",
     "the obverse/reverse setting of a rhombohedral lattice inverted"),
    ("char-mono-beta-acute", CON, "if UnitCell.bilinear(g, e1, c) > 0 { c = 0 &- c }", "if UnitCell.bilinear(g, e1, c) < 0 { c = 0 &- c }",
     "the monoclinic P cell taken with beta <= 90"),
    ("char-first-alignment-only", CON, "if score < bestScore || (score == bestScore && key.lexicographicallyPrecedes(bestKey)) {", "if bestScore == Int.max {",
     "the equivalent bases are not compared: the first is taken"),
    ("char-not-idealised", CON, "if total == .identity { cell = input } else { cell = idealised(gTotal, bravais) }", "cell = UnitCell.cell(fromMetric: gTotal)",
     "the conventional cell keeps the rounding noise of the metric"),
    ("char-automorphism-test-inverted", CON, "if abs(h - g[i][j]) > tolerance * scale[i] * scale[j] { ok = false; break search }", "if abs(h - g[i][j]) < tolerance * scale[i] * scale[j] { ok = false; break search }",
     "the lattice-symmetry test inverted"),
    # --- the structure in both settings
    ("conjugation-inverted", SET, "guard let r = transformation.conjugating(op.rotation) else { return nil }", "guard let r = transformation.conjugatedBy(op.rotation) else { return nil }",
     "R' = P R P^-1 inverted to P^-1 R P"),
    ("translation-not-transformed", SET, "translation: wrapFractional(transformation.apply(op.translation))", "translation: wrapFractional(transformation.inverseApplied(to: op.translation))",
     "t' = P t inverted to P^-1 t"),
    ("centring-images-skipped", SET, "conventional.centringTranslations.map { t in", "[SIMD3<Double>(0, 0, 0)].map { t in",
     "the lattice centring is not added to the sites"),
    ("centring-merge-too-eager", SET, "for image in images where !kept.contains(where: { matcher.matches($0.fractional, image.fractional) }) {", "for image in images where !kept.contains(where: { _ in true }) {",
     "every centring image merged into the first"),
    ("file-setting-duplicates-kept", SET, "if !fileSetting.contains(where: { $0.isEquivalent(to: back, translationTolerance: 1e-6) }) {", "if true {",
     "the centring translations are kept as operations of the file's cell"),
    ("indices-not-transformed", SET, "fileHKL.applying(transformation)", "fileHKL",
     "h' = h P not applied"),
    ("default-not-conventional", STR, "conventional: Bool = true,", "conventional: Bool = false,",
     "the default of detectingSymmetry(conventional:) false"),
    # --- tolerance and the cell list
    ("cartesian-always-matches", DET, "return cartesianLengthSquared(d) < cartesianSquared", "return true",
     "the Cartesian distance not compared"),
    ("cell-list-no-neighbour-bins", DET, "return hi - lo + 1 >= n ? (0, n) : (lo, hi - lo + 1)", "_ = hi; return (lo + Int(window[axis] * Double(n)).signum() * 0, 1)",
     "the cell list looks in one bin only (the bin too small for the window)"),
    ("cell-list-window-halved", DET, "window = SIMD3(repeating: min(0.5, t))", "window = SIMD3(repeating: min(0.5, t / 2))",
     "the cell list's search window half the tolerance"),
    ("cell-list-not-lowest-index", DET, "if j < target, !claimed[j], sites[j].z == site.z,", "if !claimed[j], sites[j].z == site.z,",
     "the target is the last candidate, not slice 1's lowest index"),
    ("cell-list-claim-ignored", DET, "if j < target, !claimed[j], sites[j].z == site.z,", "if j < target, sites[j].z == site.z,",
     "an already-claimed target may be taken twice (not a permutation)"),
]


def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def verdict_of(out):
    text = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    failed = sorted(set(re.findall(r"✘ Test (\w+)\(\)", text)))
    ran = "Test run with" in text
    if out.returncode == 0 and ran:
        return "GREEN"
    if failed:
        return "RED (%d tests: %s)" % (len(failed), ", ".join(failed))
    if ran and out.returncode != 0:
        return "RED (run failed: crash or timeout)"
    return "BUILD ERROR"


PACKAGE = '''// swift-tools-version: 6.2
import PackageDescription
let settings: [SwiftSetting] = [
    .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
    .enableUpcomingFeature("InferIsolatedConformances"),
    .enableUpcomingFeature("MemberImportVisibility"),
]
let package = Package(
    name: "DiffractoKitSlice1b",
    platforms: [.macOS("27.0")],
    products: [.library(name: "DiffractoCore", targets: ["DiffractoCore"])],
    targets: [
        .target(name: "DiffractoCore", swiftSettings: settings),
        .testTarget(name: "DiffractoCoreTests", dependencies: ["DiffractoCore"], exclude: ["Fixtures"], swiftSettings: settings),
    ],
    swiftLanguageModes: [.v6]
)
'''


def main():
    copy, build = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    wanted = set(sys.argv[3:])
    if copy.exists():
        shutil.rmtree(copy)
    (copy / "Sources").mkdir(parents=True)
    (copy / "Tests").mkdir()
    (copy / "Package.swift").write_text(PACKAGE)
    shutil.copytree(ROOT / "Sources/DiffractoCore", copy / "Sources/DiffractoCore")
    shutil.copytree(ROOT / "Tests/DiffractoCoreTests", copy / "Tests/DiffractoCoreTests")
    (copy / "Tests/DiffractoScatteringTests/Fixtures").mkdir(parents=True)
    shutil.copy(ROOT / "Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json", copy / "Tests/DiffractoScatteringTests/Fixtures/slice-2-pins.json")
    test = ["swift", "test", "-c", "release", "-j", "2", "--scratch-path", str(build)]
    base = verdict_of(run(test + ["--filter", NEW], copy))
    print("baseline (unmutated copy, the new tests): %s" % base, flush=True)
    if base != "GREEN":
        sys.exit("the unmutated copy is not green: stop")
    results = []
    for ident, rel, old, new, what in MUTANTS:
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
            out = run(test + ["--filter", NEW], copy)
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
