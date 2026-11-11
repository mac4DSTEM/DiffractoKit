//
//  RationalRelationshipTests.swift
//  The rational relationships near a measured one (DiffractoKit#81), on the published relationships of four material pairs
//  (those of TwoPassOtherMaterialsTests): what is reported when the measured relation is close, and what when it is farther off
//  than rational relationships lie from each other. The search reports; it never chooses.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoOrientation

@Suite struct RationalRelationshipTests {
    typealias Materials = TwoPassOtherMaterialsTests

    /// The search for a pair's published relationship turned by `offDeg`, and each candidate's distance to the published one.
    static func search(_ index: Int, offDeg: Double, limitDeg: Double) throws -> (found: RationalRelationshipSearch, fromPublishedDeg: [Double]) {
        let pair = Materials.pairs[index]
        let named = try OrientationRelationship.rotation(pair.first, pair.second, precipitate: pair.precipitate, matrix: pair.matrix)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(pair.matrix.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(pair.precipitate.laueClass))
        let measured = simd_double3x3(simd_quatd(angle: offDeg * .pi / 180, axis: simd_normalize(SIMD3(0.3, -0.5, 0.8)))) * named
        let found = PhaseVectorMatcher.rationalRelationships(near: measured, precipitate: pair.precipitate, matrix: pair.matrix,
                                                             matrixSymmetry: ms, precipitateSymmetry: ps, limitRad: limitDeg * .pi / 180)
        let off = found.candidates.map { PhaseVectorMatcher.relationDistance($0.relation, named, matrixSymmetry: ms, precipitateSymmetry: ps) * 180 / .pi }
        for (q, o) in zip(found.candidates, off) {
            print("RATIONAL | \(pair.name), \(offDeg) deg off, limit \(limitDeg): \(q.first.candidate) ∥ \(q.first.matrix), \(q.second.candidate) ∥ \(q.second.matrix) | index sum \(q.indexSum) | \(q.angleRad * 180 / .pi) deg from the measured, \(o) deg from the published")
        }
        return (found, off)
    }

    /// How many relationships each pair reports: (measured 0.5 degrees off, within 1 degree), (2.5 degrees off, within 4).
    /// alpha-Ti + beta-Ti 1 and 2; Ni + gamma' 1 and 1; alpha-Fe + cementite 1 and 4; Al + beta'' 1 and 4.
    static let counts: [(close: Int, loose: Int)] = [(1, 2), (1, 1), (1, 4), (1, 4)]

    /// A relation measured 0.5 degrees off its published relationship, searched within 1 degree: for every pair the published
    /// relationship is reported, alone, 0.5 degrees from the measured, and the search is not ambiguous.
    @Test(arguments: 0..<4) func aCloseMeasurementReportsThePublishedRelationshipNearest(_ index: Int) throws {
        let r = try Self.search(index, offDeg: 0.5, limitDeg: 1)
        let q = try #require(r.found.candidates.first)
        #expect(r.fromPublishedDeg[0] < 0.01, "\(Materials.pairs[index].name): the nearest is \(r.fromPublishedDeg[0]) deg from the published")
        #expect(abs(q.angleRad * 180 / .pi - 0.5) < 0.02)
        #expect(r.found.candidates.count == Self.counts[index].close, "\(Materials.pairs[index].name): \(r.found.candidates.count) candidates")
        #expect(r.found.isAmbiguous == (Self.counts[index].close > 1))
    }

    /// Measured 2.5 degrees off and searched within 4 degrees: why the search must not choose. The published relationship is
    /// among the candidates for every pair, but except for Ni + gamma' (nothing else nearby) it is not the nearest: another
    /// rational relationship, more than a degree from the published one, is nearer to the measured relation (hcp and bcc,
    /// ferrite and cementite, and a monoclinic cell in a cubic one each have several within a few degrees). The search lists
    /// them and says ambiguous.
    @Test(arguments: 0..<4) func aLooseMeasurementIsReportedAsAmbiguous(_ index: Int) throws {
        let r = try Self.search(index, offDeg: 2.5, limitDeg: 4)
        #expect(r.fromPublishedDeg.contains { $0 < 0.01 }, "the published relationship is among the candidates")
        #expect(zip(r.found.candidates, r.found.candidates.dropFirst()).allSatisfy { $0.angleRad <= $1.angleRad }, "nearest first")
        #expect(r.found.candidates.count == Self.counts[index].loose, "\(Materials.pairs[index].name): \(r.found.candidates.count) candidates")
        #expect(r.found.isAmbiguous == (index != 1))
        if index != 1 {
            #expect(try #require(r.fromPublishedDeg.first) > 1, "\(Materials.pairs[index].name): the nearest candidate is not the published relationship")
        }
    }

    /// Nothing within a limit narrower than the measurement is off; nothing for a limit that is not a positive angle below a
    /// right angle, for bounds that admit no statement, or for a matrix that is not a rotation (a scaled one, one with NaN).
    @Test func nothingWithinNothing() throws {
        #expect(try Self.search(3, offDeg: 2.5, limitDeg: 0.01).found.candidates.isEmpty)
        for limit in [Double.nan, 0, -1, 90, 200, .infinity] {
            #expect(try Self.search(3, offDeg: 2.5, limitDeg: limit).found.candidates.isEmpty, "limit \(limit)")
        }
        let pair = Materials.pairs[1]
        let named = try OrientationRelationship.rotation(pair.first, pair.second, precipitate: pair.precipitate, matrix: pair.matrix)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(pair.matrix.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(pair.precipitate.laueClass))
        func count(_ relation: simd_double3x3, maximumIndex: Int = 3, maximumIndexSum: Int = 8) -> Int {
            PhaseVectorMatcher.rationalRelationships(near: relation, precipitate: pair.precipitate, matrix: pair.matrix, matrixSymmetry: ms,
                                                     precipitateSymmetry: ps, maximumIndex: maximumIndex, maximumIndexSum: maximumIndexSum,
                                                     limitRad: 2 * .pi / 180).candidates.count
        }
        #expect(count(named) == 1, "the control: cube on cube is found for the rotation itself")
        #expect(count(2 * named) == 0, "a scaled matrix is not a rotation")
        var broken = named
        broken[0][0] = .nan
        #expect(count(broken) == 0)
        #expect(count(named, maximumIndex: 0) == 0)
        #expect(count(named, maximumIndexSum: 3) == 0)
        #expect(count(named, maximumIndex: 50) == 1, "a large index bound is cut by the index sum and returns at once")
    }
}
