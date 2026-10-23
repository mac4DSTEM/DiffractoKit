import Foundation
import Testing
@testable import DiffractoGeometry

/// Ask 6 of the apps (2026-10-22): `EllipseCalibration.FitError` (and with it `KinematicError`, `CIFImportError`, tested beside their types) is a
/// `LocalizedError` the app can show as it is.
struct FitErrorTextTests {
    @Test func everyFitErrorHasItsOwnSentence() {
        let errors: [EllipseCalibration.FitError] = [.invalidInput("a centre outside the pattern"), .insufficientSignal, .insufficientAngularCoverage(12),
                                                     .moreThanOneRing(minRadius: 20, maxRadius: 60), .didNotConverge, .invalidEllipse, .excessiveResidual(0.3), .ringWidthAtFloor]
        var seen = Set<String>()
        for e in errors {
            let text = e.localizedDescription
            #expect(text == e.errorDescription && !text.isEmpty && !text.contains("couldn’t be completed"), Comment(rawValue: text))
            seen.insert(text)
        }
        #expect(seen.count == errors.count)
    }
}
