import Testing
import DiffractoScattering

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
struct DiffractoScatteringPlaceholderTests {
    @Test func testPlannedScopeIsStated() {
        #expect(!DiffractoScattering.plannedScope.isEmpty)
    }
}
