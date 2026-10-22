import Foundation
import Testing
@testable import DiffractoCore

/// Kit asks of the apps' open items (2026-10-22) that 0.10.0 closes in DiffractoCore.
struct KitAsks0100Tests {
    /// `ReflectionFamily` has a public initialiser and keeps what it is given.
    @Test func reflectionFamilyCanBeBuiltByAConsumer() {
        let f = ReflectionFamily(g: 0.5, members: [MillerIndex(1, 1, 1), MillerIndex(-1, -1, -1)], orbits: [MillerIndex(1, 1, 1)], name: "{1 1 1}")
        #expect(f.g == 0.5 && f.d == 2 && f.multiplicity == 2 && f.name == "{1 1 1}" && f.orbits == [MillerIndex(1, 1, 1)])
        #expect(f == ReflectionFamily(g: 0.5, members: f.members, orbits: f.orbits, name: f.name))
    }

    /// `CIFImportError` is a `LocalizedError` whose text is its `description` (the kit's wording, not Foundation's "The operation couldn't be completed").
    @Test func cifImportErrorsAreLocalizedErrors() {
        let errors: [CIFImportError] = [.missingTag("_cell_length_a"), .noAtomSites, .unparsableNumber(tag: "_cell_length_a", value: "x"),
                                        .incompleteSymmetryOperationsForSymbol(symbol: "C 2/m", listed: 4, required: 8)]
        for e in errors {
            #expect(e.localizedDescription == e.description && !e.localizedDescription.contains("couldn’t be completed"), Comment(rawValue: "\(e)"))
        }
    }
}
