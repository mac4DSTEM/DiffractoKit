//
//  CIFImport+OperationChecks.swift
//  Role: the checks a CIF's symmetry-operation list must pass beyond its count, because a count cannot see wrong
//        operations: every rotation maps the cell onto itself, the list is closed (a group, modulo lattice
//        translations), and its rotations form the crystal class of the IT number the file states. Each failure is
//        refused by name (`CIFImportError.inconsistentSymmetryOperations`). pymatgen-core 2026.10.16 checks none of
//        these and expands the list as given; the kit refuses, as it refuses a short list.
//

import Foundation
import simd

/// What is wrong with a symmetry-operation list whose count passed.
public enum SymmetryListProblem: Equatable, Sendable, CustomStringConvertible {
    /// The rotation of `operation` does not map the cell onto itself (its metric changes by more than
    /// `CIFImport.Policy.operationMetricTolerance`): the list was written for another cell or another setting.
    case notACellSymmetry(operation: String)
    /// `first` applied after `second` gives `product`, which the list lacks (modulo lattice translations, compared to
    /// 1e-3; a list whose translations all lie within 0.004 of twelfths is read in twelfths): the list is not a group.
    case notClosed(first: String, second: String, product: String)
    /// The listed rotations form the crystal class `listed` (nil: none of the 32), but space group No. `spaceGroupNumber`
    /// belongs to class `required`.
    case wrongCrystalClass(spaceGroupNumber: Int, listed: String?, required: String)

    public var description: String {
        switch self {
        case .notACellSymmetry(let op):
            return "the operation \"\(op)\" does not map the cell onto itself, so the list belongs to another cell or "
                + "setting"
        case .notClosed(let first, let second, let product):
            return "\"\(first)\" after \"\(second)\" gives \"\(product)\", which the list lacks, so the list is not a "
                + "space group (an operation is wrong or missing)"
        case .wrongCrystalClass(let number, let listed, let required):
            return "its rotations form the crystal class \(listed ?? "of none of the 32 point groups"), but space group "
                + "No. \(number) belongs to class \(required)"
        }
    }
}

extension CIFImport {

    /// The crystal class (short point-group symbol) of each range of IT numbers.
    private static let crystalClasses: [(numbers: ClosedRange<Int>, symbol: String)] = [
        (1...1, "1"), (2...2, "-1"), (3...5, "2"), (6...9, "m"), (10...15, "2/m"), (16...24, "222"), (25...46, "mm2"),
        (47...74, "mmm"), (75...80, "4"), (81...82, "-4"), (83...88, "4/m"), (89...98, "422"), (99...110, "4mm"),
        (111...122, "-42m"), (123...142, "4/mmm"), (143...146, "3"), (147...148, "-3"), (149...155, "32"),
        (156...161, "3m"), (162...167, "-3m"), (168...173, "6"), (174...174, "-6"), (175...176, "6/m"),
        (177...182, "622"), (183...186, "6mm"), (187...190, "-6m2"), (191...194, "6/mmm"), (195...199, "23"),
        (200...206, "m-3"), (207...214, "432"), (215...220, "-43m"), (221...230, "m-3m"),
    ]

    /// How many elements of each rotation type a crystal class holds, in the order 1, 2, 3, 4, 6, -1, m, -3, -4, -6.
    /// The type of a rotation is read from its determinant and trace, which no change of basis alters, so the class is
    /// found in any setting. The 32 classes have 32 different signatures.
    private static let classSignatures: [[Int]: String] = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]: "1", [1, 0, 0, 0, 0, 1, 0, 0, 0, 0]: "-1",
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]: "2", [1, 0, 0, 0, 0, 0, 1, 0, 0, 0]: "m", [1, 1, 0, 0, 0, 1, 1, 0, 0, 0]: "2/m",
        [1, 3, 0, 0, 0, 0, 0, 0, 0, 0]: "222", [1, 1, 0, 0, 0, 0, 2, 0, 0, 0]: "mm2",
        [1, 3, 0, 0, 0, 1, 3, 0, 0, 0]: "mmm",
        [1, 1, 0, 2, 0, 0, 0, 0, 0, 0]: "4", [1, 1, 0, 0, 0, 0, 0, 0, 2, 0]: "-4", [1, 1, 0, 2, 0, 1, 1, 0, 2, 0]: "4/m",
        [1, 5, 0, 2, 0, 0, 0, 0, 0, 0]: "422", [1, 1, 0, 2, 0, 0, 4, 0, 0, 0]: "4mm",
        [1, 3, 0, 0, 0, 0, 2, 0, 2, 0]: "-42m", [1, 5, 0, 2, 0, 1, 5, 0, 2, 0]: "4/mmm",
        [1, 0, 2, 0, 0, 0, 0, 0, 0, 0]: "3", [1, 0, 2, 0, 0, 1, 0, 2, 0, 0]: "-3", [1, 3, 2, 0, 0, 0, 0, 0, 0, 0]: "32",
        [1, 0, 2, 0, 0, 0, 3, 0, 0, 0]: "3m", [1, 3, 2, 0, 0, 1, 3, 2, 0, 0]: "-3m",
        [1, 1, 2, 0, 2, 0, 0, 0, 0, 0]: "6", [1, 0, 2, 0, 0, 0, 1, 0, 0, 2]: "-6", [1, 1, 2, 0, 2, 1, 1, 2, 0, 2]: "6/m",
        [1, 7, 2, 0, 2, 0, 0, 0, 0, 0]: "622", [1, 1, 2, 0, 2, 0, 6, 0, 0, 0]: "6mm",
        [1, 3, 2, 0, 0, 0, 4, 0, 0, 2]: "-6m2", [1, 7, 2, 0, 2, 1, 7, 2, 0, 2]: "6/mmm",
        [1, 3, 8, 0, 0, 0, 0, 0, 0, 0]: "23", [1, 3, 8, 0, 0, 1, 3, 8, 0, 0]: "m-3",
        [1, 9, 8, 6, 0, 0, 0, 0, 0, 0]: "432", [1, 3, 8, 0, 0, 0, 6, 0, 6, 0]: "-43m",
        [1, 9, 8, 6, 0, 1, 9, 8, 6, 0]: "m-3m",
    ]

    /// The crystal class of space group `number` (1 ... 230); nil outside.
    package static func crystalClass(spaceGroup number: Int) -> String? {
        crystalClasses.first { $0.numbers.contains(number) }?.symbol
    }

    /// The crystal class the distinct rotations of `rotations` form; nil when they are not one of the 32 point groups.
    package static func crystalClass(of rotations: [IntegerMatrix]) -> String? {
        var signature = [Int](repeating: 0, count: 10)
        for rotation in Set(rotations) {
            let trace = rotation[0, 0] + rotation[1, 1] + rotation[2, 2]
            let slot: Int?
            switch (rotation.determinant, trace) {
            case (1, 3): slot = 0
            case (1, -1): slot = 1
            case (1, 0): slot = 2
            case (1, 1): slot = 3
            case (1, 2): slot = 4
            case (-1, -3): slot = 5
            case (-1, 1): slot = 6
            case (-1, 0): slot = 7
            case (-1, -1): slot = 8
            case (-1, -2): slot = 9
            default: slot = nil
            }
            guard let slot else { return nil }
            signature[slot] += 1
        }
        return classSignatures[signature]
    }

    /// The first problem of an expanding operation list, or nil: a rotation that is not an isometry of `cell` (an entry of
    /// RᵀGR − G above `metricTolerance` · sqrt(G_ii G_jj)), a product the list lacks (lists longer than
    /// `SymmetryDetection.closureCheckLimit` are not checked for closure), or a crystal class other than that of
    /// `spaceGroupNumber`.
    package static func operationListProblem(_ ops: [SymmetryOperation], cell: UnitCell, spaceGroupNumber: Int?,
                                             metricTolerance: Double) -> SymmetryListProblem? {
        // Metric: R^T G R = G for every rotation, each entry relative to sqrt(G_ii G_jj) (scale-free on every axis).
        let lattice = cell.realLattice
        let g = (0..<3).map { i in (0..<3).map { j in simd_dot(lattice[i], lattice[j]) } }
        let diagonal = [g[0][0], g[1][1], g[2][2]]
        let tolerance = metricTolerance.isFinite && metricTolerance >= 0 ? metricTolerance : 0
        if diagonal.allSatisfy({ $0.isFinite && $0 > 0 }) {
            var checked = Set<IntegerMatrix>()
            for op in ops where checked.insert(op.rotation).inserted {
                let r = op.rotation
                var worst = 0.0
                for i in 0..<3 {
                    for j in 0..<3 {
                        var sum = 0.0
                        for k in 0..<3 {
                            for l in 0..<3 where r[k, i] != 0 && r[l, j] != 0 {
                                sum += Double(r[k, i] * r[l, j]) * g[k][l]
                            }
                        }
                        worst = max(worst, abs(sum - g[i][j]) / (diagonal[i] * diagonal[j]).squareRoot())
                    }
                }
                if !(worst <= tolerance) {
                    return .notACellSymmetry(operation: op.xyzString)
                }
            }
        }

        // Closure modulo lattice translations. When every translation of the list lies within 0.004 of a twelfth (a
        // hand-rounded list: 0.33, 0.67, 0.17), all are read as those twelfths, so they compose exactly; otherwise (an
        // origin-shifted list with 0.74 or 0.11, a test shift of 0.008) they are kept and compared to 1e-3. Snapping
        // only part of a list would break the products of the rest. Rotations get small integer ids and a product
        // table; a key packs the id and the three translations in thousandths into one Int.
        func nearTwelfth(_ value: Double) -> Bool {
            guard value.isFinite else { return false }
            let scaled = (value - value.rounded(.down)) * 12
            return abs(scaled - scaled.rounded()) <= 0.004 * 12
        }
        let roundedTwelfths = ops.allSatisfy { nearTwelfth($0.translation.x) && nearTwelfth($0.translation.y)
            && nearTwelfth($0.translation.z) }
        func snapped(_ value: Double) -> Double {
            guard value.isFinite else { return 0 }
            let wrapped = value - value.rounded(.down)
            return roundedTwelfths ? (wrapped * 12).rounded() / 12 : wrapped
        }
        func milli(_ value: Double) -> Int {
            guard value.isFinite else { return 0 }
            let wrapped = value - value.rounded(.down)
            return Int((wrapped * 1000).rounded()) % 1000
        }
        var rotationIDs: [IntegerMatrix: Int] = [:]
        var rotations: [IntegerMatrix] = []
        func key(_ id: Int, _ t: SIMD3<Double>) -> Int {
            ((id * 1000 + milli(t.x)) * 1000 + milli(t.y)) * 1000 + milli(t.z)
        }
        var listed = Set<Int>()
        var distinct: [(op: SymmetryOperation, id: Int)] = []
        for raw in ops {
            let op = SymmetryOperation(rotation: raw.rotation,
                                       translation: SIMD3(snapped(raw.translation.x), snapped(raw.translation.y),
                                                          snapped(raw.translation.z)))
            let id: Int
            if let known = rotationIDs[op.rotation] { id = known } else {
                id = rotations.count
                rotationIDs[op.rotation] = id
                rotations.append(op.rotation)
            }
            if listed.insert(key(id, op.translation)).inserted { distinct.append((op, id)) }
        }
        func isListed(_ id: Int, _ t: SIMD3<Double>) -> Bool {
            if listed.contains(key(id, t)) { return true }
            let m = SIMD3(milli(t.x), milli(t.y), milli(t.z))
            for dx in -1...1 {
                for dy in -1...1 {
                    for dz in -1...1 where dx != 0 || dy != 0 || dz != 0 {
                        let packed = ((id * 1000 + (m.x + dx + 1000) % 1000) * 1000 + (m.y + dy + 1000) % 1000) * 1000
                            + (m.z + dz + 1000) % 1000
                        if listed.contains(packed) { return true }
                    }
                }
            }
            return false
        }
        func notClosed(_ first: SymmetryOperation, _ second: SymmetryOperation) -> SymmetryListProblem {
            let product = first.composed(with: second)
            return .notClosed(first: first.xyzString, second: second.xyzString,
                              product: SymmetryOperation(rotation: product.rotation,
                                                         translation: product.translationMod1).xyzString)
        }
        if distinct.count <= SymmetryDetection.closureCheckLimit {
            let productIDs = rotations.map { a in rotations.map { b in rotationIDs[a.multiplied(by: b)] } }
            for first in distinct {
                for second in distinct {
                    guard let id = productIDs[first.id][second.id] else { return notClosed(first.op, second.op) }
                    let t = first.op.rotation.apply(second.op.translation) + first.op.translation
                    if !isListed(id, t) { return notClosed(first.op, second.op) }
                }
            }
        }

        // Crystal class against the IT number.
        if let number = spaceGroupNumber, let required = crystalClass(spaceGroup: number) {
            let listedClass = crystalClass(of: rotations)
            if listedClass != required {
                return .wrongCrystalClass(spaceGroupNumber: number, listed: listedClass, required: required)
            }
        }
        return nil
    }
}
