//
//  main.swift
//  Role: uses DiffractoCore's public API from a package that depends on DiffractoKit by path (so `package` access
//        does not apply) and checks it against values a reader can verify by hand: the β″ cell of the apps, the
//        monoclinic rings, cubic gold's eight families; DiffractoGeometry's checks are in GeometryChecks.swift. Exits non-zero on any failed check. `run.sh` builds it in Debug.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry
import DiffractoOrientation

var report = Report()

report.section("UnitCell") { r in
    let cell = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)
    r.near("d(200) = a sin β / 2", cell.d(MillerIndex(2, 0, 0)), 7.31135, 6e-6)
    r.near("d(001) = c sin β", cell.d(MillerIndex(0, 0, 1)), 6.50112, 6e-6)
    r.near("d(020) = b / 2", cell.d(MillerIndex(0, 2, 0)), 2.025, 1e-12)
    r.near("volume by the closed form", cell.volume, 15.16 * 4.05 * 6.74 * sin(105.3 * .pi / 180), 1e-9)
    r.equal("metric family", cell.metricFamily(), .monoclinic(uniqueAxis: .b))
    r.equal("gLength equals 1/d", cell.gLength(MillerIndex(2, 0, -1)) * cell.d(MillerIndex(2, 0, -1)) > 0.999999999999, true)
    r.equal("cubic constructor", UnitCell.cubic(a: 4).metricFamily(), .cubic)
    r.equal("hexagonal constructor", UnitCell.hexagonal(a: 3, c: 5).metricFamily(), .hexagonal)
    r.equal("tetragonal constructor", UnitCell.tetragonal(a: 3, c: 5).metricFamily(), .tetragonal(uniqueAxis: .c))
    r.equal("tetragonal along b", UnitCell.orthorhombic(a: 3, b: 5, c: 3).metricFamily(), .tetragonal(uniqueAxis: .b))
    r.equal("tetragonal along a", UnitCell.orthorhombic(a: 5, b: 3, c: 3).metricFamily(), .tetragonal(uniqueAxis: .a))
    r.equal("orthorhombic constructor", UnitCell.orthorhombic(a: 3, b: 4, c: 5).metricFamily(), .orthorhombic)
    r.equal("rhombohedral constructor", UnitCell.rhombohedral(a: 4, alphaDeg: 70).metricFamily(), .rhombohedral)
    r.equal("triclinic constructor", UnitCell.triclinic(a: 4, b: 5, c: 6, alphaDeg: 80, betaDeg: 85, gammaDeg: 95).metricFamily(), .triclinic)
    r.equal("issues of a valid cell", cell.issues, [])
    r.equal("isValid", cell.isValid, true)
    r.near("interplanar angle", UnitCell.cubic(a: 4).interplanarAngleDeg(MillerIndex(1, 0, 0), MillerIndex(1, 1, 0)), 45, 1e-9)
    r.equal("real lattice rows", cell.realLattice.count, 3)
    r.equal("reciprocal rows", cell.reciprocalLattice.count, 3)
    r.equal("metric tensor rows", cell.metricTensor.count, 3)
    r.equal("reciprocal metric rows", cell.reciprocalMetricTensor.count, 3)
    r.near("g vector length", length(cell.gVector(MillerIndex(0, 2, 0))), 2 / 4.05, 1e-12)
    r.near("g squared from metric", cell.gSquaredFromMetric(MillerIndex(0, 2, 0)), (2 / 4.05) * (2 / 4.05), 1e-12)
}

report.section("SymmetryOperation and PointSymmetry") { r in
    let op = try SymmetryOperation(xyz: "-x,y+1/2,-z")
    r.equal("rotation rows", op.rotation, IntegerMatrix([-1, 0, 0], [0, 1, 0], [0, 0, -1]))
    r.equal("translation", op.translation, SIMD3(0, 0.5, 0))
    r.equal("xyz string", op.xyzString, "-x,y+1/2,-z")
    r.equal("apply", op.apply(SIMD3(0.1, 0.2, 0.3)), SIMD3(-0.1, 0.7, -0.3))
    r.equal("composed is identity mod 1", op.composed(with: op).isIdentity, true)
    r.equal("translation mod 1", op.composed(with: op).translationMod1, SIMD3(0, 0, 0))
    r.equal("equivalent", op.isEquivalent(to: try SymmetryOperation(xyz: "-x,y-1/2,-z")), true)
    r.equal("distinct count", SymmetryOperation.distinctCount([op, op, .identity]), 2)
    r.equal("identity constant", SymmetryOperation.identity.isIdentity, true)
    var refused = false
    do { _ = try SymmetryOperation(xyz: "x/2,y,z") } catch SymmetryOperationError.unparsable { refused = true } catch {}
    r.equal("refusal is typed", refused, true)
    let c2m = try ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z"].map { try SymmetryOperation(xyz: $0) }
    let laue = PointSymmetry.laueGroup(of: c2m)
    r.equal("Laue group order", laue.order, 4)
    r.equal("Laue class", laue.laueClass, .twoOverM)
    r.equal("crystal system", laue.laueClass?.crystalSystem, .monoclinic)
    r.equal("closed", laue.isClosed, true)
    r.equal("contains inversion", laue.containsInversion, true)
    r.equal("proper rotations", laue.properRotations.count, 2)
    r.equal("orbit of (2 0 1)", laue.orbit(of: MillerIndex(2, 0, 1)), [MillerIndex(2, 0, 1), MillerIndex(-2, 0, -1)])
    r.equal("canonical of (-2 0 1)", laue.canonical(MillerIndex(-2, 0, 1)), MillerIndex(2, 0, -1))
    r.equal("orbit name", laue.orbitName(of: MillerIndex(-2, 0, 1)), "{2 0 -1}")
    r.equal("holohedry of cubic", PointSymmetry.holohedry(of: .cubic).order, 48)
    r.equal("hexagonal name", PointSymmetry.holohedry(of: .hexagonal).orbitName(of: MillerIndex(-1, -1, 0)), "{1 1 -2 0}")
    r.equal("closure of generators", PointSymmetry(closureOf: [.inversion]).order, 2)
    r.equal("rotations as given", PointSymmetry(rotations: [.identity, .inversion]).contains(.inversion), true)
    r.equal("Laue class symbol", LaueClass.mThreeBarM.symbol, "m-3m")
    r.equal("Laue class order", LaueClass.sixOverMmm.order, 24)
    r.equal("Laue class Miller-Bravais", LaueClass.threeBar.usesMillerBravais, true)
    r.equal("Laue class from point group", LaueClass(pointGroupSymbol: "32"), .threeBarM)
    r.equal("MillerIndex applying", MillerIndex(1, 2, 3).applying(.inversion), MillerIndex(-1, -2, -3))
    r.equal("MillerIndex order", MillerIndex(1, 0, 0) > MillerIndex(0, 1, 0), true)
    r.equal("MillerIndex i", MillerIndex(1, 1, 0).millerBravaisI, -2)
    r.equal("IntegerMatrix order", IntegerMatrix.inversion.order, 2)
    r.equal("IntegerMatrix determinant", IntegerMatrix.identity.determinant, 1)
    r.equal("Element symbol", Element.symbol(forAtomicNumber: 12), "Mg")
    r.equal("Element number", Element.atomicNumber(forSymbol: "Si"), 14)
    r.equal("Element from CIF", Element.normalizedSymbol(fromCIF: "FE3+"), "Fe")
    r.equal("Element count", Element.symbols.count, Element.maximumAtomicNumber)
}

report.section("CrystalStructure, CIFImport, ReflectionList") { r in
    let cif = """
    data_gold
    _cell_length_a 4.0782(3)
    _cell_length_b 4.0782
    _cell_length_c 4.0782
    _cell_angle_alpha 90
    _cell_angle_beta 90
    _cell_angle_gamma 90
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    Au1 Au 0 0 0
    Au2 Au 0.5 0.5 0
    Au3 Au 0.5 0 0.5
    Au4 Au 0 0.5 0.5
    """
    let doc = try CIFImport.read(cif, fileBaseName: "gold")
    r.equal("display name", doc.displayName, "gold")
    r.equal("data block", doc.dataBlockName, "gold")
    r.equal("declared symbol", doc.declaredSpaceGroupSymbol, nil)
    r.equal("declared number", doc.declaredSpaceGroupNumber, nil)
    r.equal("half-step", doc.coarsestHalfStep, 0)
    r.equal("had operations", doc.hadSymmetryOperations, false)
    let s = doc.structure
    r.equal("sites", s.sites.count, 4)
    r.equal("site label", s.sites[0].label, "Au1")
    r.equal("site Z", s.sites[0].z, 79)
    r.equal("site occupancy", s.sites[0].occupancy, 1)
    r.equal("formula", s.formula, "Au4")
    r.equal("symmetry source", s.symmetrySource, .detected)
    r.equal("operations", s.operations.count, 192)
    r.equal("Laue class", s.laueClass, .mThreeBarM)
    r.equal("crystal system", s.crystalSystem, .cubic)
    r.equal("Laue group order", s.laueGroup.order, 48)
    r.equal("space group number", s.spaceGroupNumber, nil)
    r.equal("name", s.name, "gold")
    r.equal("unresolved", s.unresolvedRotations.isEmpty, true)
    r.equal("usable", s.isUsable, true)
    r.equal("validation issues", s.validationIssues.count, 0)
    r.equal("close contact", s.shortestCloseContact, nil)
    r.equal("close contact limit", CrystalStructure.closeContactLimit, 0.5)
    let list = try s.reflections(kMax: 1.2)
    r.equal("k_max", list.kMax, 1.2)
    r.equal("allowed (100)", ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: s.operations), false)
    r.equal("first reflection is {111}", list.allowed.first?.hkl.h.magnitude, 1)
    r.near("first d", list.allowed.first?.d ?? 0, 4.0782 / 3.0.squareRoot(), 1e-12)
    r.equal("reflection fields", list.reflections.first.map { $0.isAllowed || !$0.isAllowed && length($0.g) == $0.gLength } ?? false, true)
    let families = list.families()
    r.equal("families", families.map(\.name), ["{1 1 1}", "{2 0 0}", "{2 2 0}", "{3 1 1}", "{2 2 2}", "{4 0 0}", "{3 3 1}", "{4 2 0}"])
    r.equal("multiplicity", families[0].multiplicity, 8)
    r.equal("members", families[0].members.count, 8)
    r.equal("orbits", families[0].orbits, [MillerIndex(1, 1, 1)])
    r.near("family d", families[0].d, 4.0782 / 3.0.squareRoot(), 1e-12)
    r.near("family g", families[0].g, 3.0.squareRoot() / 4.0782, 1e-12)
    r.equal("lattice families are more", list.families(allowedOnly: false).count > families.count, true)
    r.equal("list Laue class", list.laueClass, .mThreeBarM)
    r.equal("constants", (ReflectionList.maximumKMax, ReflectionList.familyTolerance) == (4, 1e-9), true)
    var refused: ReflectionListError?
    do { _ = try s.reflections(kMax: 0) } catch let e as ReflectionListError { refused = e } catch {}
    r.equal("k_max refusal", refused, .invalidKMax(0))
    var missing: CIFImportError?
    do { _ = try CIFImport.read("data_x\n_cell_length_a 1\n") } catch let e as CIFImportError { missing = e } catch {}
    r.equal("CIF refusal is typed", missing, .missingTag("_cell_length_b"))
    r.equal("space group order", CIFImport.spaceGroupOrder(225), 192)
    r.equal("required count", CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 105, gammaDeg: 90), 8)
    r.equal("half-step rule", CIFImport.coordinateHalfStep("0.333"), 1e-3)
    let built = CrystalStructure(cell: .cubic(a: 4), sites: [AtomSite(z: 26, fractional: [0, 0, 0])],
                                 operations: [.identity, SymmetryOperation(rotation: .inversion)], spaceGroupNumber: 2, name: "x")
    r.equal("given operations", built.symmetrySource, .operationsGiven)
    r.equal("given Laue class", built.laueClass, .oneBar)
    let detected = CrystalStructure.detectingSymmetry(cell: .cubic(a: 2.8665), sites: [
        AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.5])])
    r.equal("detected operations", detected.operations.count, 96)
    let detection = SymmetryDetection.detect(sites: detected.sites, cell: detected.cell)
    r.equal("detection result", detection.operations.count, 96)
    r.equal("detection candidates", detection.candidates.order, 48)
    r.equal("detection tolerance", SymmetryDetection.positionTolerance, 5e-3)
    r.equal("fractional helpers", fractionalAxisDelta(0.99, 0.01) < 0.021 && wrapFractional(SIMD3(-0.25, 1.5, 0)) == SIMD3(0.75, 0.5, 0), true)
}

report.section("DiffractoScattering (slice 2)") { r in
    let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: fcc.map { AtomSite(z: 79, fractional: $0) })
    let list = try gold.kinematicReflections(kMax: 1.2)
    r.equal("Au rings after the cut", list.families.count, 8)
    r.equal("Au first ring", list.families.first?.family.name, "{1 1 1}")
    r.equal("cut default", StructureFactors.defaultCut, 1e-4)
    r.equal("kept count", list.kept.count, 112)
    r.equal("all >= kept", list.all.count >= list.kept.count, true)
    r.equal("k_max and cut", list.kMax == 1.2 && list.cut == 1e-4, true)
    r.equal("no operation-absent kept", list.keptButOperationAbsent.isEmpty, true)
    let f111 = try unwrap(StructureFactors.compute(MillerIndex(1, 1, 1), in: gold))
    let fe = try unwrap(ElectronScatteringFactors.f(z: 79, gSquared: f111.gLength * f111.gLength))
    r.near("F(111) = 4 f_e / V", f111.real, 4 * fe / gold.cell.volume, 1e-12)
    r.near("|F| and |F|^2", f111.magnitude * f111.magnitude, f111.intensity, 1e-15)
    r.near("d", f111.d, 4.0782 / 3.0.squareRoot(), 1e-12)
    r.equal("imaginary part zero", abs(f111.imaginary) < 1e-15, true)
    r.equal("g vector", f111.g == gold.cell.gVector(MillerIndex(1, 1, 1)), true)
    r.equal("ring intensity positive", (list.families.first?.intensity ?? 0) > 0, true)
    r.equal("reflection in kinematic list", list.kept.first?.reflection.hkl.h.magnitude, 1)
    r.equal("table size", ElectronScatteringFactors.parameters.count, ElectronScatteringFactors.maximumAtomicNumber)
    r.equal("supported", ElectronScatteringFactors.isSupported(z: 103) && !ElectronScatteringFactors.isSupported(z: 104), true)
    let custom: FormFactor = { _, _ in 1 }
    r.equal("custom form factor", try gold.kinematicReflections(kMax: 1.2, cut: 0, formFactor: custom).kept.isEmpty, false)
    var refused: KinematicError?
    do { _ = try gold.kinematicReflections(kMax: 1.2, cut: -1) } catch let e as KinematicError { refused = e } catch {}
    r.equal("cut refusal typed", refused, .invalidCut(-1))
}

runGeometryChecks(&report)

exit(report.finish())

