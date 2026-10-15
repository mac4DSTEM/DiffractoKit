//
//  KinematicErrorTextTests.swift
//  Registration 0.0.4 section 5 (tier 2): `KinematicError` is a `LocalizedError` whose `localizedDescription` is its
//  `description`.
//

import Foundation
import Testing
import DiffractoCore
import DiffractoScattering

struct KinematicErrorTextTests {

    @Test func testLocalizedDescriptionIsTheDescription() {
        let errors: [KinematicError] = [
            .unsupportedElements([104, 118]),
            .reflectionList(.invalidKMax(5)),
            .invalidCut(-1),
            .invalidStructure(["atomic sites too close"]),
            .nonFiniteStructureFactor(MillerIndex(1, 1, 1)),
        ]
        for error in errors {
            #expect(!error.description.isEmpty)
            #expect(error.localizedDescription == error.description, "\(error)")
            #expect((error as Error).localizedDescription == error.description, "\(error) through Error")
        }
        // The nested list error keeps its own text.
        #expect(KinematicError.reflectionList(.invalidKMax(5)).localizedDescription
                == ReflectionListError.invalidKMax(5).localizedDescription)
    }

    @Test func testAThrownErrorShowsTheKitsText() throws {
        let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: [AtomSite(z: 79, fractional: [0, 0, 0])])
        let error = try #require(throws: KinematicError.self) { try gold.kinematicReflections(kMax: 1.2, cut: -1) }
        #expect(error == .invalidCut(-1))
        #expect(error.localizedDescription.hasPrefix("The |F| cut must be"))
    }
}
