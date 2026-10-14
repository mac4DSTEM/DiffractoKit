import XCTest
import DiffractoGeometry

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
final class DiffractoGeometryPlaceholderTests: XCTestCase {
    func testPlannedScopeIsStated() {
        XCTAssertFalse(DiffractoGeometry.plannedScope.isEmpty)
    }
}
