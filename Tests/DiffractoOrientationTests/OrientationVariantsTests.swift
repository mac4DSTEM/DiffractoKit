import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// beta'' (Mg5Si6, C2/m) in Al: the 12 orientation-relationship variants against the planted truth of tools/beta-variants/planted.py (Andersen et al.
/// 1998: b || [001]Al, a || [3 2 0]Al, c in the (001) plane). Closed forms, no Python reference: validation "none".
struct OrientationVariantsTests {
    static let matrixGroup = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
    static let precipitateGroup = CrystalOrientationSymmetry(laueClass: .twoOverM)

    /// The planted R0: columns a/|a| = [3 2 0]/sqrt(13), b = [0 0 1], their cross product.
    static var plantedR0: simd_double3x3 {
        let ea = SIMD3(3.0, 2.0, 0.0) / 13.0.squareRoot(), eb = SIMD3(0.0, 0.0, 1.0)
        return simd_double3x3(columns: (ea, eb, simd_cross(ea, eb)))
    }

    static func relation() throws -> simd_double3x3 {
        try OrientationRelationship.rotation(
            OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(3, 2, 0))),
            OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
            precipitate: CrystalPresets.betaDoublePrime, matrix: CrystalPresets.aluminum)
    }

    static func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var g = 0.0
        for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
        return g
    }
    static func angle(_ m: simd_double3x3) -> Double { let q = simd_quatd(m); return 2 * atan2(simd_length(q.imag), abs(q.real)) }

    @Test func theStatedRelationBuildsThePlantedRotation() throws {
        #expect(Self.gap(try Self.relation(), Self.plantedR0) < 1e-12)
        // c of the precipitate lies 105.25 degrees from a: [1 -3 0] is the stated direction of c (planted.py asserts the angle)
        let r = try Self.relation()
        let cell = CrystalPresets.betaDoublePrime.cell
        let c = r * simd_normalize(cell.realLattice[2])
        let target = simd_normalize(SIMD3(1.0, -3.0, 0.0))
        #expect(acos(simd_dot(c, target)) * 180 / .pi < 0.1)
    }

    @Test func thereAreTwelveVariantsFourPerNeedleAxis() throws {
        let set = try OrientationVariantSet(matrix: Self.plantedR0, matrixSymmetry: Self.matrixGroup, precipitateSymmetry: Self.precipitateGroup)
        #expect(set.variants.count == 12)
        // b of the precipitate (column 1) in the Al frame is the needle axis: +-[100], +-[010], +-[001], four variants for each cubic axis
        var perAxis = [0, 0, 0]
        for v in set.variants {
            let b = v.columns.1
            let axis = (0 ..< 3).max { abs(b[$0]) < abs(b[$1]) }!
            #expect(abs(abs(b[axis]) - 1) < 1e-12)
            perAxis[axis] += 1
        }
        #expect(perAxis == [4, 4, 4])
        // the variants are 12 different classes, and every equivalent of the relation belongs to one of them
        for (i, a) in set.variants.enumerated() {
            for (j, b) in set.variants.enumerated() where i < j {
                #expect(!Self.precipitateGroup.operators.contains { Self.gap(a * $0, b) < 1e-9 })
            }
        }
        let all = Self.matrixGroup.operators.map { $0 * Self.plantedR0 }
        #expect(all.allSatisfy { r in set.variants.contains { v in Self.precipitateGroup.operators.contains { Self.gap(v * $0, r) < 1e-9 } } })
    }

    /// Variant-to-variant disorientation (2/m on both sides): only 0, 90, 120 and 180 degrees, 60 / 48 / 24 ordered pairs; the relation's own
    /// reduced disorientation (432 on the left, 2/m on the right) is arctan(2/3).
    @Test func variantDisorientationsAndTheRelationsOwnAngle() throws {
        let set = try OrientationVariantSet(matrix: Self.plantedR0, matrixSymmetry: Self.matrixGroup, precipitateSymmetry: Self.precipitateGroup)
        var counts: [Int: Int] = [:]
        for (i, a) in set.variants.enumerated() {
            for (j, b) in set.variants.enumerated() where i != j {
                let angle = Self.precipitateGroup.operators.flatMap { s1 in
                    Self.precipitateGroup.operators.map { s2 in Self.angle(s1 * a.transpose * b * s2) }
                }.min()! * 180 / .pi
                let rounded = Int(angle.rounded())
                #expect(abs(angle - Double(rounded)) < 1e-9 && [0, 90, 120, 180].contains(rounded), Comment(rawValue: "\(angle)"))
                counts[rounded, default: 0] += 1
            }
        }
        #expect(counts[90] == 60 && counts[120] == 48 && counts[180] == 24 && counts[0] == nil)
        let reduced = try Misorientation(Self.plantedR0, source: Self.precipitateGroup, target: Self.matrixGroup).reduce()
        #expect(abs(reduced.angleRad - atan2(2.0, 3.0)) < 1e-12)
    }

    @Test func nearestVariantFindsThePlantedVariantOfNoisyOrientations() throws {
        let set = try OrientationVariantSet(matrix: Self.plantedR0, matrixSymmetry: Self.matrixGroup, precipitateSymmetry: Self.precipitateGroup)
        var state: UInt64 = 35
        func uniform() -> Double {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return Double(state >> 11) / Double(1 << 53) * 2 - 1
        }
        for _ in 0 ..< 48 {
            let k = Int((uniform() * 0.5 + 0.5) * 12) % 12
            let axis = simd_normalize(SIMD3(uniform(), uniform(), uniform()))
            let delta = simd_double3x3(simd_quatd(angle: 0.5 * .pi / 180, axis: axis))
            let s = Self.precipitateGroup.operators[uniform() > 0 ? 0 : 1]
            let measured = set.variants[k] * s * delta
            let hit = try #require(set.nearestVariant(of: measured))
            #expect(hit.index == k)
            #expect(abs(hit.disorientationRad * 180 / .pi - 0.5) < 1e-9 + 0.0001 && (hit.nextRad ?? 0) * 180 / .pi > 89)
        }
        #expect(set.nearestVariant(of: simd_double3x3(rows: [SIMD3(1, 1, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)])) == nil)
    }

    /// ACOM / phase-map orientations: matrix orientation M_m and precipitate orientation M_p (columns = lab axes in crystal coordinates).
    @Test func observedRelationFromTwoOrientationMatrices() throws {
        let set = try OrientationVariantSet(matrix: Self.plantedR0, matrixSymmetry: Self.matrixGroup, precipitateSymmetry: Self.precipitateGroup)
        let lab = simd_double3x3(simd_quatd(angle: 1.1, axis: simd_normalize(SIMD3(1, 2, 3))))     // matrix orientation: lab to Al crystal
        let variant = set.variants[7]
        let precipitate = variant.transpose * lab                                                     // lab to precipitate crystal, so M_m M_p^T = variant
        let hit = try #require(set.nearestVariant(matrixOrientation: lab, precipitateOrientation: precipitate))
        #expect(hit.index == 7 && hit.disorientationRad < 1e-9)
    }

    @Test func refusals() throws {
        let notRotation = simd_double3x3(rows: [SIMD3(2, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)])
        #expect(throws: OrientationSymmetryError.invalidRotation) { try OrientationVariantSet(matrix: notRotation, matrixSymmetry: Self.matrixGroup, precipitateSymmetry: Self.precipitateGroup) }
        #expect(throws: OrientationSymmetryError.invalidRotation) { try Misorientation(notRotation, source: Self.matrixGroup, target: Self.matrixGroup) }
        #expect(throws: OrientationSymmetryError.invalidRotation) { try Orientation(notRotation, symmetry: Self.matrixGroup) }
        // two statements about the same axis do not fix the rotation
        #expect(throws: OrientationSymmetryError.invalidDirection) {
            try OrientationRelationship.rotation(
                OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
                OrientationRelationship(candidate: .direction(SIMD3(0, 2, 0)), matrix: .direction(SIMD3(0, 0, 1))),
                precipitate: CrystalPresets.betaDoublePrime, matrix: CrystalPresets.aluminum)
        }
    }
}
