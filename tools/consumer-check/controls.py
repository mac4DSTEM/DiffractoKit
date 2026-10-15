#!/usr/bin/env python3
"""consumer-check negative controls (run through `run.sh --controls`; Python 3 standard library only).

SpectroKit's tools/consumer-check/controls.py pattern. The consumer package (tools/consumer-check) compiles only while
every symbol it uses is `public`. This proves it: for each such symbol ONE declaration is reverted from `public` to
`package` in a COPY of the repository (rsync into <root>/.build/consumer-check/mutants/<name>, never the working tree), the
consumer is built against the copy (its Package.swift takes the dependency by path; ONE warm build directory serves every
control, so each control is an incremental build), and the build must FAIL with a compiler diagnostic that names that
symbol. A control is caught only when
  * the unedited copy builds (the baseline, so a failure is the edit's), and
  * the edit applied exactly once, and
  * the build exits non-zero, and the compiler's diagnostic is the expected one:
      "access":  the consumer is refused, `error: '<member>' is inaccessible due to 'package' protection level`, and the
                 note beside it points at the edited declaration (`DiffractoX.<Owner>.<member>: note: '<member>(labels)'
                 declared here`), so two declarations with one name are told apart;
      "decl":    a whole type becomes `package`: the consumer is refused (`'<Type>' is inaccessible due to 'package'
                 protection level`, `cannot find type '<Type>' in scope`, `'<Type>' is not a member type of ...`), or, when
                 a public signature of the library itself mentions the type, the library no longer compiles
                 (`... cannot be declared public because ... uses a package type`, the code frame under it naming the type);
      "used":    a declaration whose refusal has no stable name in the diagnostic (an initialiser, an overload): the
                 consumer's build must fail with an error in one of the consumer's own source files matching a pattern;
      "library": a declaration the library itself refuses before any consumer compiles (it witnesses a public protocol
                 requirement, or a public default argument uses it): an error in the edited library file matching a pattern.
  * Declarations that cannot be reverted from `public` to `package` are listed in SKIPPED with the reason (an enum case, a
    protocol conformance): they are not counted.

usage: controls.py <repository root> [--dry] [<symbol> ...]     (no symbol: all of them)
       --dry: no build; check that every edit target occurs exactly once in its file, and print the counts.
environment: CONSUMER_CHECK_JOBS=N passes `-j N` to every swift build.
exit:  0 only when every control is caught.
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

ARGS = [a for a in sys.argv[1:] if a != "--dry"]
DRY = "--dry" in sys.argv[1:]
ROOT = os.path.abspath(ARGS[0])
ONLY = ARGS[1:]
NAME = os.path.basename(ROOT)
MUTANTS = os.path.join(ROOT, ".build", "consumer-check", "mutants")
COPY = os.path.join(MUTANTS, NAME)
LOGS = os.path.join(MUTANTS, "logs")
CONSUMER = os.path.join(COPY, "tools", "consumer-check")
JOBS = os.environ.get("CONSUMER_CHECK_JOBS", "")

CORE = "Sources/DiffractoCore/"
SCAT = "Sources/DiffractoScattering/"
GEO = "Sources/DiffractoGeometry/"
ANSI = re.compile(r"\x1b\[[0-9;]*m")
WITNESS = r"must be declared public because it matches a requirement in public protocol"
DEFAULT_ARGUMENT = r"static property 'default' is package and cannot be referenced from a default argument value|is package and cannot be referenced from a default argument value"


def module_of(rel):
    return re.match(r"Sources/(\w+)/", rel).group(1)


def access(rel, owner, member, old):
    """The declaration `old` (its first `public`) becomes `package`; the consumer is refused at the member's use."""
    return dict(kind="access", symbol=f"{owner}.{member}", rel=rel, old=old, owner=owner, member=member,
                module=module_of(rel))


def decl(rel, symbol, old, word=None):
    """A whole type becomes `package`. `word` is the name the diagnostic carries (the last component by default)."""
    return dict(kind="decl", symbol=symbol, rel=rel, old=old, word=word or symbol.split(".")[-1])


def used(rel, symbol, old, pattern="inaccessible due to 'package' protection level|no exact matches|cannot convert|"
                                    "incorrect argument|cannot be constructed|missing argument"):
    """No stable name in the diagnostic (an initialiser, an overload): an error in the consumer's own sources matching
    `pattern`."""
    return dict(kind="used", symbol=symbol, rel=rel, old=old, pattern=pattern)


def library(rel, symbol, old, pattern, file=None):
    """The library itself refuses the edit: an error in the file `file` (the edited file `rel` by default) matching `pattern`."""
    return dict(kind="library", symbol=symbol, rel=rel, old=old, pattern=pattern, file=file or os.path.basename(rel))


UC = CORE + "UnitCell.swift"
CS = CORE + "CrystalStructure.swift"
CIF = CORE + "CIFImport.swift"
RL = CORE + "ReflectionList.swift"
LC = CORE + "LaueClass.swift"
PS = CORE + "PointSymmetry.swift"
MI = CORE + "MillerIndex.swift"
EL = CORE + "Element.swift"
SD = CORE + "SymmetryDetection.swift"
SOP = CORE + "SymmetryOperation.swift"
EPOL = GEO + "EllipseCalibration+Policy.swift"
LR = CORE + "LatticeReduction.swift"
CC = CORE + "ConventionalCell.swift"
CSET = CORE + "ConventionalSetting.swift"
SF = SCAT + "StructureFactors.swift"
ESF = SCAT + "ElectronScatteringFactors.swift"
EW = GEO + "ElectronWavelength.swift"
RC = GEO + "RingCalibration.swift"
RCF = GEO + "RingCalibration+Fit.swift"
RCQ = GEO + "RingCalibration+Q.swift"
RE = GEO + "RingEvidence.swift"
RI = GEO + "RingIndexing.swift"
ZA = GEO + "ZoneAxis.swift"
DP = GEO + "DiffractionPattern.swift"
EC = GEO + "EllipseCalibration.swift"
ET = GEO + "EllipseTransform.swift"

# Each `old` occurs exactly once in its file (checked at run time, `--dry` checks without a build).
CONTROLS = [
    # ---- 0.0.4's new public symbols (docs/registration-0.0.4.md) -----------------------------------------------------
    # RingCalibration.fitRing / calibrating / q / radius
    access(RCF, "RingCalibration", "fitRing", "public static func fitRing(pattern:"),
    access(RCF, "RingCalibration", "calibrating", "public static func calibrating("),
    access(RCQ, "RingCalibration", "q", "public static func q(radius:"),
    # `RingCalibration.radius` is also an instance property: with the static function package the call resolves to it.
    used(RCQ, "RingCalibration.radius(q:cameraLengthPixels:wavelengthAngstroms:)", "public static func radius(q:",
         r"instance member 'radius' cannot be used on type 'RingCalibration'"),
    # RingRadiusFit: the type and its fields, init, radius, isUsable
    decl(RCF, "RingRadiusFit", "public struct RingRadiusFit: Sendable, Equatable {"),
    access(RCF, "RingRadiusFit", "fit", "public let fit: EllipseCalibrationFit?\n    public let report: RingCalibration.FitReport?\n"
           "    public let refusals: [RingCalibration.Refusal]\n\n    public init"),
    access(RCF, "RingRadiusFit", "report", "public let report: RingCalibration.FitReport?\n"
           "    public let refusals: [RingCalibration.Refusal]\n\n    public init"),
    access(RCF, "RingRadiusFit", "refusals", "public let refusals: [RingCalibration.Refusal]\n\n    public init"),
    used(RCF, "RingRadiusFit.init", "public init(fit: EllipseCalibrationFit?, report: RingCalibration.FitReport?, refusals: [RingCalibration.Refusal]) {"),
    access(RCF, "RingRadiusFit", "radius", "public var radius: Double? { report?.radius }"),
    access(RCF, "RingRadiusFit", "isUsable", "public var isUsable: Bool"),
    # The four new refusal cases (enum cases) are SKIPPED; their name and reason are existing members; errorDescription is
    # the LocalizedError witness.
    # FitReport's new init parameters and fields
    access(RC, "RingCalibration.FitReport", "annulusInnerRadius", "public var annulusInnerRadius: Double?"),
    access(RC, "RingCalibration.FitReport", "annulusOuterRadius", "public var annulusOuterRadius: Double?"),
    access(RC, "RingCalibration.FitReport", "ringEvidence", "public var ringEvidence: RingEvidence?"),
    used(RC, "RingCalibration.FitReport.init(annulusInnerRadius:annulusOuterRadius:ringEvidence:)",
         "public init(model: EllipseCalibrationModel, semiAxisA: Double, semiAxisB: Double, occupiedAngularBins: Int,"),
    # RingEvidence: the type, every field, the init, the bounds, the predicates
    decl(RE, "RingEvidence", "public struct RingEvidence: Sendable, Equatable {"),
    access(RE, "RingEvidence", "amplitude", "public let amplitude: Double"),
    access(RE, "RingEvidence", "secondaryPeak", "public let secondaryPeak: Double?"),
    access(RE, "RingEvidence", "maximumSecondaryPeak", "public static let maximumSecondaryPeak = 0.25"),
    access(RE, "RingEvidence", "isSingleRing", "public var isSingleRing: Bool"),
    access(RE, "RingEvidence", "sectorAmplitudes", "public let sectorAmplitudes: [Double?]"),
    access(RE, "RingEvidence", "sectorAmplitudeErrors", "public let sectorAmplitudeErrors: [Double?]"),
    access(RE, "RingEvidence", "sectorShifts", "public let sectorShifts: [Double?]"),
    access(RE, "RingEvidence", "occupiedBins", "public let occupiedBins: Int"),
    access(RE, "RingEvidence", "innerDrop", "public let innerDrop: Double?"),
    access(RE, "RingEvidence", "outerDrop", "public let outerDrop: Double?"),
    access(RE, "RingEvidence", "radiusShift", "public let radiusShift: Double?"),
    access(RE, "RingEvidence", "centreShift", "public let centreShift: Double?"),
    access(RE, "RingEvidence", "centreShiftError", "public let centreShiftError: Double?"),
    access(RE, "RingEvidence", "ellipseShift", "public let ellipseShift: Double?"),
    access(RE, "RingEvidence", "ellipseShiftError", "public let ellipseShiftError: Double?"),
    used(RE, "RingEvidence.init", "public init(amplitude: Double, sectorAmplitudes:"),
    access(RE, "RingEvidence", "sectorCount", "public static let sectorCount = 36"),
    access(RE, "RingEvidence", "minimumSectorPixels", "public static let minimumSectorPixels = 8"),
    access(RE, "RingEvidence", "minimumAmplitudeFraction", "public static let minimumAmplitudeFraction = 0.5"),
    access(RE, "RingEvidence", "minimumSignificance", "public static let minimumSignificance = 3.0"),
    access(RE, "RingEvidence", "minimumOccupiedBins", "public static let minimumOccupiedBins = 30"),
    access(RE, "RingEvidence", "minimumDrop", "public static let minimumDrop = 0.25"),
    access(RE, "RingEvidence", "maximumShift", "public static let maximumShift = 0.25"),
    access(RE, "RingEvidence", "shiftSignificance", "public static let shiftSignificance = 4.0"),
    access(RE, "RingEvidence", "isPeak", "public var isPeak: Bool"),
    access(RE, "RingEvidence", "geometryHolds", "public var geometryHolds: Bool"),
    access(RE, "RingEvidence", "passes", "public var passes: Bool"),
    # EllipseCalibration.ringEvidence and EllipseCalibrationFit.ringEvidence
    access(RE, "EllipseCalibration", "ringEvidence", "public static func ringEvidence(pattern:"),
    access(EC, "EllipseCalibrationFit", "ringEvidence", "public let ringEvidence: RingEvidence?"),
    used(EC, "EllipseCalibrationFit.init(ringEvidence:)", "public init(centerQX: Double, centerQY: Double, a: Double, b: Double, theta: Double,"),
    # requireRingEvidence / progress / isCancelled sit on these declarations: reverting one removes the whole call
    access(EC, "EllipseCalibration", "fitBestAvailable", "public static func fitBestAvailable("),
    access(EC, "EllipseCalibration", "fitAmorphousRing", "public static func fitAmorphousRing("),
    access(EC, "EllipseCalibration", "fit1D", "public static func fit1D("),
    access(RCF, "RingCalibration", "fit", "public static func fit(pattern:"),
    # UnitCell.MetricFamily (the type of .tetragonal(uniqueAxis:)); LocalizedError's errorDescription on the three types
    decl(UC, "UnitCell.MetricFamily", "public enum MetricFamily: Sendable, Equatable {"),
    library(CIF, "CIFImportError.errorDescription", "public var errorDescription: String? { description }",
            r"property 'errorDescription' " + WITNESS + r" 'LocalizedError'"),
    library(RL, "ReflectionListError.errorDescription", "public var errorDescription: String? { description }",
            r"property 'errorDescription' " + WITNESS + r" 'LocalizedError'"),
    library(SF, "KinematicError.errorDescription", "public var errorDescription: String? { description }",
            r"property 'errorDescription' " + WITNESS + r" 'LocalizedError'"),
    library(RC, "RingCalibration.Refusal.errorDescription", "public var errorDescription: String? { reason }",
            r"property 'errorDescription' " + WITNESS + r" 'LocalizedError'"),
    access(RC, "RingCalibration.Refusal", "name", "public var name: String {"),
    access(RC, "RingCalibration.Refusal", "reason", "public var reason: String {"),
    access(RC, "RingCalibration.Refusal", "withholdsDeltaQ", "public var withholdsDeltaQ: Bool {"),

    # ---- DiffractoCore ----------------------------------------------------------------------------------------------
    decl(UC, "UnitCell", "public struct UnitCell: Sendable, Equatable {"),
    access(UC, "UnitCell", "monoclinic", "public static func monoclinic(a:"),
    access(UC, "UnitCell", "d", "public func d(_ hkl"),
    access(UC, "UnitCell", "metricFamily", "public func metricFamily("),
    access(UC, "UnitCell", "gVector", "public func gVector("),
    decl(CS, "CrystalStructure", "public struct CrystalStructure: Sendable, Equatable {"),
    access(CS, "CrystalStructure", "detectingSymmetry", "public static func detectingSymmetry("),
    access(CS, "CrystalStructure", "laueClass", "public let laueClass: LaueClass?"),
    access(CS, "CrystalStructure", "operations", "public let operations: [SymmetryOperation]"),
    access(CS, "CrystalStructure", "formula", "public var formula: String {"),
    access(CS, "CrystalStructure", "reflections", "public func reflections(kMax: Double)"),
    decl(CS, "AtomSite", "public struct AtomSite: Sendable, Equatable {"),
    access(CS, "AtomSite", "z", "public var z: Int"),
    access(CS, "AtomSite", "label", "public var label: String?"),
    used(CS, "AtomSite.init", "public init(z: Int, fractional: SIMD3<Double>, occupancy: Double = 1, label: String? = nil) {"),
    decl(CIF, "CIFImport", "public enum CIFImport {"),
    access(CIF, "CIFImport", "read", "public static func read(_ cifText"),
    access(CIF, "CIFImport", "spaceGroupOrder", "public static func spaceGroupOrder("),
    decl(CIF, "CIFDocument", "public struct CIFDocument: Sendable, Equatable {"),
    access(CIF, "CIFDocument", "structure", "public let structure: CrystalStructure"),
    access(CIF, "CIFDocument", "displayName", "public let displayName: String"),
    decl(RL, "ReflectionList", "public struct ReflectionList: Sendable, Equatable {"),
    access(RL, "ReflectionList", "allowed", "public var allowed: [Reflection]"),
    access(RL, "ReflectionList", "isAllowed", "public static func isAllowed("),
    access(RL, "ReflectionList", "maximumKMax", "public static let maximumKMax"),
    decl(LC, "LaueClass", "public enum LaueClass: String, Sendable, Equatable, CaseIterable, CustomStringConvertible {"),
    access(LC, "LaueClass", "symbol", "public var symbol: String"),
    access(LC, "LaueClass", "order", "public var order: Int {"),
    used(LC, "LaueClass.init(pointGroupSymbol:)", "public init?(pointGroupSymbol raw: String) {"),
    decl(PS, "PointSymmetry", "public struct PointSymmetry: Sendable, Equatable {"),
    access(PS, "PointSymmetry", "holohedry", "public static func holohedry(of family"),
    access(PS, "PointSymmetry", "order", "public var order: Int"),
    access(PS, "PointSymmetry", "orbit", "public func orbit(of hkl: MillerIndex)"),
    decl(MI, "MillerIndex", "public struct MillerIndex: Sendable, Hashable, Comparable, CustomStringConvertible {"),
    access(MI, "MillerIndex", "h", "public let h: Int"),
    access(MI, "MillerIndex", "applying", "public func applying("),
    decl(EL, "Element", "public enum Element {"),
    access(EL, "Element", "symbol", "public static func symbol(forAtomicNumber"),
    access(EL, "Element", "normalizedSymbol", "public static func normalizedSymbol(fromCIF"),
    # a public default argument of CrystalStructure.detectingSymmetry uses SymmetryDetection: the library refuses the edit
    library(SD, "SymmetryDetection", "public enum SymmetryDetection {",
            r"enum 'SymmetryDetection' is package and cannot be referenced from a default argument value", "CrystalStructure.swift"),
    access(SD, "SymmetryDetection", "positionTolerance", "public static let positionTolerance"),
    access(SD, "SymmetryDetection", "detect", "public static func detect("),

    # ---- slice 1b's new public symbols (docs/registration-slice-1b.md) -------------------------------------------------
    # Niggli reduction: ReducedCell, UnitCell.niggliReduced / niggliRelativeEpsilon
    decl(LR, "ReducedCell", "public struct ReducedCell: Sendable, Equatable {"),
    access(LR, "ReducedCell", "cell", "public let cell: UnitCell\n    /// P with determinant +1"),
    access(LR, "ReducedCell", "transformation", "public let transformation: IntegerMatrix\n    /// ε in Å²"),
    access(LR, "ReducedCell", "epsilon", "public let epsilon: Double"),
    access(LR, "ReducedCell", "converged", "public let converged: Bool"),
    access(LR, "ReducedCell", "steps", "public let steps: Int"),
    access(LR, "UnitCell", "niggliReduced", "public func niggliReduced("),
    # the default argument of niggliReduced(relativeEpsilon:) uses the constant: the library refuses the edit
    library(LR, "UnitCell.niggliRelativeEpsilon", "public static let niggliRelativeEpsilon",
            r"static property 'niggliRelativeEpsilon' is package and cannot be referenced from a default argument value"),
    # The conventional cell: BravaisLattice, ConventionalCell, UnitCell.conventionalCell / conventionalCellTolerance
    decl(CC, "BravaisLattice", "public enum BravaisLattice: String, Sendable, Equatable, CaseIterable, CustomStringConvertible {"),
    library(CC, "BravaisLattice.description", "public var description: String { rawValue }",
            r"property 'description' " + WITNESS + r" 'CustomStringConvertible'"),
    access(CC, "BravaisLattice", "centring", "public var centring: Character"),
    access(CC, "BravaisLattice", "holohedryOrder", "public var holohedryOrder: Int"),
    access(CC, "BravaisLattice", "centringMultiplicity", "public var centringMultiplicity: Int {\n        switch centring"),
    decl(CC, "ConventionalCell", "public struct ConventionalCell: Sendable, Equatable {"),
    access(CC, "ConventionalCell", "cell", "public let cell: UnitCell\n    public let bravais"),
    access(CC, "ConventionalCell", "bravais", "public let bravais: BravaisLattice"),
    access(CC, "ConventionalCell", "transformation", "public let transformation: IntegerMatrix\n    public let reduced"),
    access(CC, "ConventionalCell", "reduced", "public let reduced: ReducedCell"),
    access(CC, "ConventionalCell", "centringTranslations", "public let centringTranslations: [SIMD3<Double>]"),
    access(CC, "ConventionalCell", "uniqueAxis", "public let uniqueAxis: UnitCell.Axis?"),
    access(CC, "ConventionalCell", "setting", "public let setting: String\n\n    public var centringMultiplicity"),
    access(CC, "ConventionalCell", "centringMultiplicity", "public var centringMultiplicity: Int { transformation.determinant }"),
    access(CC, "UnitCell", "conventionalCell", "public func conventionalCell(tolerance"),
    # the default argument of conventionalCell(tolerance:) uses the constant
    library(CC, "UnitCell.conventionalCellTolerance", "public static let conventionalCellTolerance",
            r"static property 'conventionalCellTolerance' is package and cannot be referenced from a default argument value"),
    # ConventionalSetting: the type, its fields, the shortcuts, the structure and the four conversions
    decl(CSET, "ConventionalSetting", "public struct ConventionalSetting: Sendable, Equatable {"),
    access(CSET, "ConventionalSetting", "conventionalCell", "public let conventionalCell: ConventionalCell"),
    access(CSET, "ConventionalSetting", "sites", "public let sites: [AtomSite]"),
    access(CSET, "ConventionalSetting", "operations", "public let operations: [SymmetryOperation]"),
    access(CSET, "ConventionalSetting", "laueGroup", "public let laueGroup: PointSymmetry"),
    access(CSET, "ConventionalSetting", "laueClass", "public let laueClass: LaueClass?"),
    access(CSET, "ConventionalSetting", "translationUncertainty", "public let translationUncertainty: Double"),
    access(CSET, "ConventionalSetting", "cell", "public var cell: UnitCell { conventionalCell.cell }"),
    access(CSET, "ConventionalSetting", "bravais", "public var bravais: BravaisLattice { conventionalCell.bravais }"),
    access(CSET, "ConventionalSetting", "transformation", "public var transformation: IntegerMatrix { conventionalCell.transformation }"),
    access(CSET, "ConventionalSetting", "structure", "public var structure: CrystalStructure {"),
    access(CSET, "ConventionalSetting", "conventionalIndices", "public func conventionalIndices(of"),
    access(CSET, "ConventionalSetting", "fileIndices", "public func fileIndices(ofConventional"),
    access(CSET, "ConventionalSetting", "fileOperation", "public func fileOperation(of"),
    access(CSET, "ConventionalSetting", "conventionalOperation", "public func conventionalOperation(of"),
    # CrystalStructure.Setting, .setting, .conventionalSetting
    decl(CS, "CrystalStructure.Setting", "public enum Setting: Sendable, Equatable {"),
    access(CS, "CrystalStructure", "setting", "public let setting: Setting"),
    access(CS, "CrystalStructure", "conventionalSetting", "public let conventionalSetting: ConventionalSetting?"),
    # SymmetryDetection.Tolerance and .Search (public default arguments of detectingSymmetry and detect use both)
    decl(SD, "SymmetryDetection.Tolerance", "public enum Tolerance: Sendable, Equatable {"),
    library(SD, "SymmetryDetection.Tolerance.default", "public static let `default`: Tolerance",
            r"static property 'default' is package and cannot be referenced from a default argument value", "CrystalStructure.swift"),
    access(SD, "SymmetryDetection.Tolerance", "proposedCartesian", "public static let proposedCartesian: Tolerance"),
    decl(SD, "SymmetryDetection.Search", "public enum Search: Sendable, Equatable {"),

    # ---- DiffractoScattering ----------------------------------------------------------------------------------------
    # StructureFactors.swift's public default arguments use ElectronScatteringFactors.f: the library refuses the edits
    library(ESF, "ElectronScatteringFactors", "public enum ElectronScatteringFactors {",
            r"enum 'ElectronScatteringFactors' is package and cannot be referenced from a default argument value", "StructureFactors.swift"),
    library(ESF, "ElectronScatteringFactors.f", "public static func f(z:",
            r"static method 'f\(z:gSquared:\)' is package and cannot be referenced from a default argument value", "StructureFactors.swift"),
    access(ESF, "ElectronScatteringFactors", "isSupported", "public static func isSupported("),
    decl(SF, "KinematicReflectionList", "public struct KinematicReflectionList: Sendable, Equatable {"),
    access(SF, "KinematicReflectionList", "kept", "public let kept: [KinematicReflection]"),
    access(SF, "KinematicReflectionList", "families", "public let families: [KinematicFamily]"),
    decl(SF, "KinematicReflection", "public struct KinematicReflection: Sendable, Equatable {"),
    access(SF, "KinematicReflection", "structureFactor", "public let structureFactor: StructureFactor"),
    decl(SF, "KinematicFamily", "public struct KinematicFamily: Sendable, Equatable {"),
    access(SF, "KinematicFamily", "family", "public let family: ReflectionFamily\n    /// Σ|F|²"),
    access(SF, "KinematicFamily", "intensity", "public let intensity: Double\n}\n\n/// Every reflection"),

    # ---- DiffractoGeometry ------------------------------------------------------------------------------------------
    decl(EW, "ElectronWavelength", "public enum ElectronWavelength {"),
    access(EW, "ElectronWavelength", "angstroms", "public static func angstroms(volts"),
    access(EW, "ElectronWavelength", "planckConstant", "public static let planckConstant"),
    decl(RC, "RingCalibration", "public struct RingCalibration: Sendable, Equatable {"),
    access(RC, "RingCalibration", "calibrate", "public static func calibrate("),
    access(RC, "RingCalibration", "refusals", "public static func refusals("),
    # `twoTheta` is also an instance property: with the static function package the call resolves to it.
    used(RC, "RingCalibration.twoTheta(wavelengthAngstroms:d:)", "public static func twoTheta(",
         r"instance member 'twoTheta' cannot be used on type 'RingCalibration'"),
    access(RC, "RingCalibration", "deltaQ", "public let deltaQ: Double"),
    decl(RC, "RingCalibration.FitReport", "public struct FitReport: Sendable, Equatable {"),
    decl(RC, "RingCalibration.Refusal", "public enum Refusal: Error, Equatable, Sendable, LocalizedError {"),
    decl(RCF, "RingFitCalibration", "public struct RingFitCalibration: Sendable, Equatable {"),
    access(RCF, "RingFitCalibration", "fit", "public let fit: EllipseCalibrationFit?\n    public let report: RingCalibration.FitReport?\n"
           "    public let refusals: [RingCalibration.Refusal]\n    public let calibration"),
    access(RCF, "RingFitCalibration", "calibration", "public let calibration: RingCalibration?"),
    decl(RI, "RingIndexing", "public enum RingIndexing {"),
    library(RI, "RingIndexing.defaultCoverage", "public static let defaultCoverage",
            r"static property 'defaultCoverage' is package and cannot be referenced from a default argument value"),
    access(RI, "RingIndexing", "reach", "public static func reach("),
    used(RI, "RingIndexing.index(in: KinematicReflectionList)", "public static func index(q: Double, sigma: Double, in list: KinematicReflectionList,"),
    used(RI, "RingIndexing.index(in: [KinematicFamily])", "public static func index(q: Double, sigma: Double, in families: [KinematicFamily], kMax: Double,"),
    used(RI, "RingIndexing.index(in: [ReflectionFamily])", "public static func index(q: Double, sigma: Double, in families: [ReflectionFamily], kMax: Double,"),
    decl(ZA, "ZoneAxis", "public struct ZoneAxis: Sendable, Hashable, CustomStringConvertible {"),
    access(ZA, "ZoneAxis", "lowIndex", "public static let lowIndex"),
    access(ZA, "ZoneAxis", "cartesian", "public func cartesian(in cell"),
    access(ZA, "ZoneAxis", "detectorBasis", "public static func detectorBasis("),
    access(ZA, "ZoneAxis", "project", "public static func project("),
    decl(DP, "DiffractionPattern", "public struct DiffractionPattern: Sendable, Equatable {"),
    access(DP, "DiffractionPattern", "width", "public let width: Int"),
    used(DP, "DiffractionPattern.init", "public init(width: Int, height: Int, pixels: [Float]) {"),
    decl(EC, "EllipseCalibration", "public enum EllipseCalibration {"),
    decl(EC, "EllipseCalibration.FitError", "public enum FitError: LocalizedError, Equatable {"),
    decl(EC, "EllipseCalibrationFit", "public struct EllipseCalibrationFit: Sendable, Equatable {"),
    access(EC, "EllipseCalibrationFit", "a", "public let a: Double"),
    access(EC, "EllipseCalibrationFit", "profile", "public let profile: EllipseProfileParameters?"),
    decl(ET, "EllipseTransform", "public struct EllipseTransform: Sendable, Equatable {"),
    access(ET, "EllipseTransform", "accepted", "public static func accepted("),
    access(ET, "EllipseTransform", "uncorrectedOffset", "public func uncorrectedOffset("),
    used(ET, "EllipseTransform.init", "public init?(a: Double?, b: Double?, theta: Double?) {"),

    # ---- DiffractoKit 0.2.0: the Policy types (docs/registration-0.2.0.md), retainingAll, the tolerance domain -----------
    # A `static let default` that a public initialiser uses as a default argument: the library refuses the edit.
    # ReflectionList.Policy
    decl(RL, "ReflectionList.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(RL, "ReflectionList.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(RL, "ReflectionList.Policy", "maximumKMax", "public var maximumKMax: Double"),
    access(RL, "ReflectionList.Policy", "maximumReflections", "public var maximumReflections: Int"),
    access(RL, "ReflectionList.Policy", "familyTolerance", "public var familyTolerance: Double"),
    used(RL, "ReflectionList.Policy.init", "public init(maximumKMax: Double = ReflectionList.maximumKMax,"),
    access(RL, "ReflectionList", "policy", "public let policy: Policy"),
    # CrystalStructure.ValidationPolicy
    decl(CS, "CrystalStructure.ValidationPolicy", "public struct ValidationPolicy: Sendable, Equatable {", "ValidationPolicy"),
    library(CS, "CrystalStructure.ValidationPolicy.default", "public static let `default` = ValidationPolicy()", DEFAULT_ARGUMENT),
    access(CS, "CrystalStructure.ValidationPolicy", "closeContactLimit", "public var closeContactLimit: Double"),
    access(CS, "CrystalStructure", "isUsable", "public func isUsable(under policy"),
    access(CS, "CrystalStructure", "shortestCloseContact", "public func shortestCloseContact(using policy"),
    access(CS, "CrystalStructure", "validationIssues", "public func validationIssues(using policy"),
    # SymmetryDetection.Policy, the tolerance's domain and refusal, SymmetryOperation(xyz:maximumMagnitude:)
    decl(SD, "SymmetryDetection.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(SD, "SymmetryDetection.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(SD, "SymmetryDetection.Policy", "occupancyTolerance", "public var occupancyTolerance: Double"),
    access(SD, "SymmetryDetection.Tolerance", "refusal", "public func refusal(for cell: UnitCell)"),
    access(SD, "SymmetryDetection.Tolerance", "isValid", "public func isValid(for cell: UnitCell)"),
    access(SD, "SymmetryDetection.Tolerance", "minimumWindow", "public static let minimumWindow"),
    access(SD, "SymmetryDetection.Result", "refusal", "public let refusal: String?"),
    access(CS, "CrystalStructure", "symmetryRefusal", "public let symmetryRefusal: String?"),
    used(SOP, "SymmetryOperation.init(xyz:maximumMagnitude:)", "public init(xyz raw: String, maximumMagnitude"),
    library(SOP, "SymmetryOperation.defaultMaximumMagnitude", "public static let defaultMaximumMagnitude: Double", DEFAULT_ARGUMENT),
    # CIFImport.Policy
    decl(CIF, "CIFImport.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(CIF, "CIFImport.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(CIF, "CIFImport.Policy", "positionTolerance", "public var positionTolerance: SymmetryDetection.Tolerance"),
    access(CIF, "CIFImport.Policy", "maximumMergeDistance", "public var maximumMergeDistance: Double"),
    # EllipseCalibration.Policy
    decl(EPOL, "EllipseCalibration.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(EPOL, "EllipseCalibration.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(EPOL, "EllipseCalibration.Policy", "angularBinCount", "public var angularBinCount: Int"),
    access(EPOL, "EllipseCalibration.Policy", "degeneracyBoundBins", "public var degeneracyBoundBins: Int"),
    access(EPOL, "EllipseCalibration.Policy", "isValid", "public var isValid: Bool"),
    access(EPOL, "EllipseCalibration.Policy", "evidence", "public var evidence: RingEvidence.Policy"),
    used(EPOL, "EllipseCalibration.Policy.init", "public init(\n            angularBinCount: Int = 36,"),
    # RingEvidence.Policy and its judgements
    decl(RE, "RingEvidence.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(RE, "RingEvidence.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(RE, "RingEvidence.Policy", "minimumDrop", "public var minimumDrop: Double"),
    access(RE, "RingEvidence.Policy", "isValid", "public var isValid: Bool"),
    used(RE, "RingEvidence.isPeak(_:)", "public func isPeak(_ policy: Policy) -> Bool", "inaccessible due to 'package' protection level|cannot be used|no exact matches"),
    used(RE, "RingEvidence.passes(_:)", "public func passes(_ policy: Policy) -> Bool", "inaccessible due to 'package' protection level|cannot be used|no exact matches"),
    # RingCalibration.Policy
    decl(RC, "RingCalibration.Policy", "public struct Policy: Sendable, Equatable {", "Policy"),
    library(RC, "RingCalibration.Policy.default", "public static let `default` = Policy()", DEFAULT_ARGUMENT),
    access(RC, "RingCalibration.Policy", "minimumOccupiedAngularBins", "public var minimumOccupiedAngularBins: Int"),
    access(RC, "RingCalibration.Policy", "maximumResidual", "public var maximumResidual: Double"),
    used(RC, "RingCalibration.Refusal.reason(policy:)", "public func reason(policy: RingCalibration.Policy)", "inaccessible due to 'package' protection level|no exact matches|cannot be used"),
    # RingIndexing's minimum margin; kinematicReflections (retainingAll:, policy:, validation:)
    library(RI, "RingIndexing.minimumMargin", "public static let minimumMargin: Double", DEFAULT_ARGUMENT),
    used(SF, "CrystalStructure.kinematicReflections(retainingAll:)", "public func kinematicReflections(kMax: Double, cut: Double"),
]

# What cannot be reverted from `public` to `package`, and why (printed, not counted).
SKIPPED = [
    ("RingCalibration.Refusal.ringOutsideAnnulus", "an enum case takes its enum's access: it cannot be declared package on its own"),
    ("RingCalibration.Refusal.notARingPeak", "an enum case takes its enum's access (the enum itself is controlled above)"),
    ("RingCalibration.Refusal.ringGeometryInconsistent", "an enum case takes its enum's access (the enum itself is controlled above)"),
    ("RingCalibration.Refusal.secondRingInAnnulus", "an enum case takes its enum's access (the enum itself is controlled above)"),
    ("UnitCell.MetricFamily.tetragonal(uniqueAxis:)", "an enum case takes its enum's access (the enum itself is controlled above)"),
    ("BravaisLattice.aP ... cF", "enum cases take their enum's access (the enum itself is controlled above)"),
    ("CrystalStructure.Setting.file / .conventional", "enum cases take their enum's access (the enum itself is controlled above)"),
    ("SymmetryDetection.Tolerance.fractional / .cartesian", "enum cases take their enum's access (the enum itself is controlled above)"),
    ("SymmetryDetection.Search.cellList / .bruteForce", "enum cases take their enum's access (the enum itself is controlled above)"),
    ("CrystalStructure.detectingSymmetry(conventional:tolerance:search:)", "parameters have no access of their own; the function is controlled above (detectingSymmetry)"),
    ("SymmetryDetection.detect(tolerance:search:)", "parameters have no access of their own; the function is controlled above (detect)"),
    ("CIFImportError: LocalizedError", "a protocol conformance has no access of its own; its witness errorDescription is controlled above"),
    ("ReflectionListError: LocalizedError", "a protocol conformance has no access of its own; its witness errorDescription is controlled above"),
    ("KinematicError: LocalizedError", "a protocol conformance has no access of its own; its witness errorDescription is controlled above"),
]


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def build(log_name):
    """swift build -c debug of the consumer in the copy: (exit status, the output with colour codes removed)."""
    command = ["swift", "build", "-c", "debug"] + (["-j", JOBS] if JOBS else [])
    result = subprocess.run(command, cwd=CONSUMER, capture_output=True, text=True, timeout=1800)
    text = ANSI.sub("", result.stdout + result.stderr)
    with open(os.path.join(LOGS, log_name + ".log"), "w", encoding="utf-8") as f:
        f.write(text)
    return result.returncode, text


def short(line):
    """A diagnostic line without the long path before the file name."""
    return re.sub(r"^\S*/([^/\s]+\.swift)", r"\1", line.strip())


def judge(control, status, text):
    """(caught, the evidence line, why not)."""
    if status == 0:
        return False, None, "the build SUCCEEDED with the declaration reverted to package"
    kind = control["kind"]
    if kind == "access":
        member, owner, module = control["member"], control["owner"], control["module"]
        error = re.search(rf"^(\S+\.swift:\d+:\d+): error: '{re.escape(member)}' is inaccessible due to 'package' protection level$",
                          text, re.M)
        declared = re.search(rf"^{re.escape(module)}\.{re.escape(owner)}\.{re.escape(member)}\S*: note: "
                             rf"'{re.escape(member)}(?:\([^']*\))?' declared here$", text, re.M)
        if not error:
            return False, None, f"no `error: '{member}' is inaccessible due to 'package' protection level` line"
        if not declared:
            return False, None, (f"no note `{module}.{owner}.{member}...: note: '{member}...' declared here` "
                                 f"(a different declaration?)")
        return True, short(error.group(0)), None
    if kind == "decl":
        word = control["word"]
        refused = re.search(rf"^(\S+\.swift:\d+:\d+): error: ('(?:\w+\.)*{re.escape(word)}' is inaccessible due to 'package' protection level"
                            rf"|cannot find (?:type )?'(?:\w+\.)*{re.escape(word)}' in scope"
                            rf"|'{re.escape(word)}' is not a member type of .*)$", text, re.M)
        if refused:
            return True, short(refused.group(0)), None
        # else the library's own public signature mentions the type: the code frame under the error names it
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if "uses a package type" in line and "cannot be declared public because" in line:
                frame = "\n".join(lines[max(0, i - 6):i + 8])
                if re.search(rf"\b{re.escape(word)}\b", frame):
                    return True, f"library: {short(line)[:120]}  [frame names {word}]", None
        return False, None, f"neither the consumer refused `{word}` nor a public signature of the library naming it fails"
    if kind == "library":
        for m in re.finditer(r"^(\S+\.swift):\d+:\d+: error: (.*)$", text, re.M):
            if m.group(1).endswith("/" + control["file"]) and re.search(control["pattern"], m.group(2)):
                return True, f"library: {short(m.group(0))[:170]}", None
        return False, None, f"no library error in {control['file']} matching /{control['pattern']}/"
    # kind == "used"
    for m in re.finditer(r"^(\S+/tools/consumer-check/Sources/\S+\.swift:\d+:\d+): error: (.*)$", text, re.M):
        if re.search(control["pattern"], m.group(2)):
            return True, short(m.group(0))[:200], None
    return False, None, f"no error in the consumer's sources matching /{control['pattern']}/"


def main():
    selected = [c for c in CONTROLS if not ONLY or c["symbol"] in ONLY]
    unknown = [s for s in ONLY if s not in {c["symbol"] for c in CONTROLS}]
    if unknown:
        print("consumer-check: unknown control(s): " + ", ".join(unknown), file=sys.stderr)
        return 2
    symbols = [c["symbol"] for c in CONTROLS]
    if len(set(symbols)) != len(symbols):
        print("consumer-check: duplicate control(s): " + ", ".join(sorted({s for s in symbols if symbols.count(s) > 1})), file=sys.stderr)
        return 2
    started = time.time()

    if DRY:
        bad = 0
        for c in selected:
            with open(os.path.join(ROOT, c["rel"]), encoding="utf-8") as f:
                n = f.read().count(c["old"])
            if n != 1:
                bad += 1
                print(f"  FAIL  {c['symbol']}: the edit target occurs {n} times in {c['rel']}: {c['old']!r}")
        print(f"consumer-check: dry run: {len(selected) - bad} of {len(selected)} edit targets occur exactly once")
        return 1 if bad else 0

    # The repository's own files are never written; their hashes are compared at the end.
    originals = {rel: sha256(os.path.join(ROOT, rel)) for rel in sorted({c["rel"] for c in CONTROLS})}

    # A fresh copy: a directory under <root>/.build only, then the sources without .build, .git and References.
    assert MUTANTS.startswith(os.path.join(ROOT, ".build") + os.sep) and NAME, MUTANTS
    shutil.rmtree(MUTANTS, ignore_errors=True)
    os.makedirs(COPY)
    os.makedirs(LOGS)
    subprocess.run(["rsync", "-a", "--exclude", ".build", "--exclude", ".git", "--exclude", "References",
                    ROOT + "/", COPY + "/"], check=True)

    status, text = build("baseline")
    if status != 0:
        print("consumer-check: the UNEDITED copy does not build, so no control means anything:", file=sys.stderr)
        print("\n".join(text.splitlines()[-25:]), file=sys.stderr)
        return 2
    print(f"consumer-check: controls build in {COPY} (baseline: the unedited copy builds, one warm build directory)")

    caught = 0
    for number, control in enumerate(selected, 1):
        symbol, rel = control["symbol"], control["rel"]
        path = os.path.join(COPY, rel)
        with open(os.path.join(ROOT, rel), encoding="utf-8") as f:
            pristine = f.read()
        if pristine.count(control["old"]) != 1:
            print(f"  FAIL  {symbol}: the edit target occurs {pristine.count(control['old'])} times in {rel}, not once"
                  f" (re-aim it): {control['old']!r}")
            continue
        edited = pristine.replace(control["old"], control["old"].replace("public", "package", 1))
        with open(path, "w", encoding="utf-8") as f:
            f.write(edited)
        try:
            status, text = build(f"{number:03d}-{re.sub(r'[^A-Za-z0-9]+', '-', symbol)}")
        finally:
            with open(path, "w", encoding="utf-8") as f:
                f.write(pristine)
        ok, evidence, why = judge(control, status, text)
        if ok:
            caught += 1
            print(f"  ok    {symbol}: {evidence}", flush=True)
        else:
            print(f"  FAIL  {symbol}: {why} (build exit {status}; log {LOGS})", flush=True)

    for symbol, reason in SKIPPED:
        print(f"  skip  {symbol}: {reason}")

    # The copy is restored and the repository untouched.
    for rel, digest in originals.items():
        if sha256(os.path.join(ROOT, rel)) != digest:
            print(f"  FAIL  the repository's {rel} changed during the controls")
            return 1
        if sha256(os.path.join(COPY, rel)) != digest:
            print(f"  FAIL  the copy's {rel} was not restored")
            return 1
    total = len(selected)
    print(f"consumer-check: negative controls: {caught} of {total} caught ({time.time() - started:.0f} s; "
          f"{len(SKIPPED)} declarations skipped as not revertible; the repository's sources untouched, the copy restored)")
    return 0 if caught == total else 1


if __name__ == "__main__":
    sys.exit(main())
