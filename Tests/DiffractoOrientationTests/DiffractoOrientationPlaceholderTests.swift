import Testing
import DiffractoOrientation

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
struct DiffractoOrientationPlaceholderTests {
    @Test func testPlannedScopeIsStated() {
        #expect(!DiffractoOrientation.plannedScope.isEmpty)
    }
}
