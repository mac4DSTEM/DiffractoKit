import XCTest
import DiffractoScattering

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
final class DiffractoScatteringPlaceholderTests: XCTestCase {
    func testPlannedScopeIsStated() {
        XCTAssertFalse(DiffractoScattering.plannedScope.isEmpty)
    }
}
