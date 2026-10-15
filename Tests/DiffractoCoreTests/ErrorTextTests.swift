//
//  ErrorTextTests.swift
//  Registration 0.0.4 section 5 (tier 2): the typed errors of DiffractoCore are `LocalizedError`s whose
//  `localizedDescription` is their `description`, so an app that shows `error.localizedDescription` shows the kit's text
//  rather than Foundation's "The operation couldn't be completed. (DiffractoCore.CIFImportError error 0.)".
//

import Foundation
import Testing
import DiffractoCore

struct ErrorTextTests {

    @Test func testCIFImportErrorLocalizedDescriptionIsItsDescription() {
        let errors: [CIFImportError] = [
            .missingTag("_cell_length_a"),
            .unparsableNumber(tag: "_cell_length_b", value: "four"),
            .noAtomSites,
            .truncatedLoop(firstTag: "_atom_site_label", columns: 5, values: 7),
            .invalidStructure(["atoms 0.3 Å apart"]),
        ]
        for error in errors {
            #expect(!error.description.isEmpty)
            #expect(error.localizedDescription == error.description, "\(error)")
            #expect((error as Error).localizedDescription == error.description, "\(error) through Error")
        }
    }

    @Test func testReflectionListErrorLocalizedDescriptionIsItsDescription() {
        let errors: [ReflectionListError] = [
            .invalidKMax(5), .invalidKMax(.nan), .cellTooLargeOrInvalid(estimatedReflections: 3e9),
            .cellTooLargeOrInvalid(estimatedReflections: .infinity),
        ]
        for error in errors {
            #expect(!error.description.isEmpty)
            #expect(error.localizedDescription == error.description, "\(error)")
            #expect((error as Error).localizedDescription == error.description, "\(error) through Error")
        }
    }

    /// The text is thrown by the library itself, not only built by hand.
    @Test func testAThrownErrorShowsTheKitsText() throws {
        let error = try #require(throws: ReflectionListError.self) {
            try ReflectionList(cell: .cubic(a: 4), operations: [.identity], laueGroup: PointSymmetry(closureOf: []),
                               laueClass: nil, kMax: -1)
        }
        #expect(error.localizedDescription.hasPrefix("k_max must be a number"))
        #expect(!error.localizedDescription.contains("couldn’t be completed"))
    }
}
