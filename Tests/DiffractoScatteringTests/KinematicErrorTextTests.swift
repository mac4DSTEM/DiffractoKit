import Foundation
import Testing
@testable import DiffractoScattering

struct KinematicErrorTextTests {
    @Test func theLocalizedTextIsTheKitsDescription() {
        let errors: [KinematicError] = [.unsupportedElements([118]), .invalidCut(-1), .invalidStructure(["a close contact"])]
        for e in errors { #expect(e.localizedDescription == e.description && !e.localizedDescription.contains("couldn’t be completed")) }
    }
}
