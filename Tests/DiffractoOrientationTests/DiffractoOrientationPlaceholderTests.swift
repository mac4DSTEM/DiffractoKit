import XCTest
import DiffractoOrientation

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
final class DiffractoOrientationPlaceholderTests: XCTestCase {
    func testPlannedScopeIsStated() {
        XCTAssertFalse(DiffractoOrientation.plannedScope.isEmpty)
    }
}
