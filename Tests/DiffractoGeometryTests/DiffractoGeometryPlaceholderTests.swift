import Testing
import DiffractoGeometry

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
struct DiffractoGeometryPlaceholderTests {
    @Test func testPlannedScopeIsStated() {
        #expect(!DiffractoGeometry.plannedScope.isEmpty)
    }
}
