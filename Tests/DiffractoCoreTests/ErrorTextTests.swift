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

    /// Registration 0.2.0, ask 2 (tier 3, a text test): a cell too large for the default limit says how many reflections it
    /// holds, what the default limit is, and which policy field sets another; an invalid cell says it is not a lattice.
    @Test func testTheTooLargeRefusalNamesItsNumbersAndTheField() throws {
        let text = ReflectionListError.cellTooLargeOrInvalid(estimatedReflections: 3e9).description
        #expect(text.contains("3e+09"), "\(text)")
        #expect(text.contains("1e7"), "\(text)")
        #expect(text.contains("ReflectionList.Policy.maximumReflections"), "\(text)")
        #expect(text.contains("maximumEnumeration"), "\(text)")
        #expect(ReflectionList.maximumReflections == 10_000_000)
        #expect(ReflectionListError.cellTooLargeOrInvalid(estimatedReflections: .infinity).description.contains("not a valid finite lattice"))
        // thrown by the library itself, for a 60 A cell at the default k_max cap
        let error = try #require(throws: ReflectionListError.self) {
            try ReflectionList(cell: .cubic(a: 60), operations: [.identity], laueGroup: PointSymmetry(closureOf: []),
                               laueClass: nil, kMax: 4)
        }
        #expect(error.localizedDescription.contains("ReflectionList.Policy.maximumReflections"))
    }
}
