import Foundation
import Testing
@testable import DiffractoGeometry

/// `LinearScale.assign` with kinematically forbidden shells that double diffraction reaches (mac4DSTEM#274): they may match a later ring,
/// never the first, and they don't push a crystal's own shells out of the window a hypothesis uses. Planted closed forms, plus one real
/// scan when its radii are on this machine.
@Suite("LinearScale.assign, double-diffraction shells")
struct DoubleDiffractionShellTests {
    static let a = 5.4309
    /// Diamond Si's allowed shells by h² + k² + l²: 111, 220, 311, 400, 331, 422, 511/333, 440, 531, 620.
    static let si: [Double] = [3, 8, 11, 16, 19, 24, 27, 32, 35, 40].map { Double($0).squareRoot() / a }
    static let siLabels = ["{111}", "{220}", "{311}", "{400}", "{331}", "{422}", "{511}", "{440}", "{531}", "{620}"]
    /// The forbidden shells g1 + g2 reaches: 200, 222, 420, 600/442, 622.
    static let siDoubleDiffraction: [Double] = [4, 12, 20, 36, 44].map { Double($0).squareRoot() / a }
    static let deltaQ = 0.0205

    static func positions(_ rings: [Double]) -> [[Double]] {
        LinearScaleAssignmentTests.positions(rings.map { $0 / deltaQ }, count: 40, noise: 0.05)
    }

    static func assign(_ radii: [[Double]], doubleDiffraction: [Double]) throws -> LinearScaleAssignment? {
        try LinearScale.assign(radiiPerPosition: radii, shellLengths: si, labels: siLabels, crystalName: "Si",
                               doubleDiffractionShells: doubleDiffraction)
    }

    /// Si [110]: {111}, {002} (double diffraction) and {220}, the {002} ring 1 % above its kinematic radius as measured on the real scan
    /// (18.16 px pooled against 17.98 from {111}): its ratio to {111} is then that of {311} to {220}, and without the forbidden shell the
    /// sequence is read from {220} or refused. The offset is what makes this test discriminate: at the exact kinematic radius the kept
    /// shells alone already give {111}. Mutation: ignore `doubleDiffractionShells` - red.
    @Test func diamond110IsReadFromItsFirstShell() throws {
        let radii = Self.positions([Self.si[0], 1.01 * Self.siDoubleDiffraction[0], Self.si[1]])
        let without: LinearScaleAssignment?? = try? Self.assign(radii, doubleDiffraction: [])
        #expect(without.flatMap { $0 }?.label != "{111}", "without the forbidden shell the sequence is read from another shell (or refused)")
        let found = try #require(try Self.assign(radii, doubleDiffraction: Self.siDoubleDiffraction))
        #expect(found.label == "{111}")
        #expect(found.matchedRings == 3)
        #expect(abs(found.referenceRadiusInvAngstrom / found.observedRadiiPixels[0] / Self.deltaQ - 1) < 0.01)
        #expect(found.nextShellInvAngstrom == Self.siDoubleDiffraction[0], "the second ring is the {002} it matched")
        #expect(found.rmsMisfit < 0.012)
    }

    /// Si [001]: {220}, {400}, {440}. A forbidden shell must not start a reading ({200} would fit 1 : √2 : 2 too) nor push {440} out of
    /// the window. Mutation: forbidden shells merged into the shell list - {200} first, red.
    @Test func diamond001KeepsItsAssignment() throws {
        let radii = Self.positions([Self.si[1], Self.si[3], Self.si[7]])
        let without = try #require(try Self.assign(radii, doubleDiffraction: []))
        let with = try #require(try Self.assign(radii, doubleDiffraction: Self.siDoubleDiffraction))
        #expect(without.label == "{220}")
        #expect(with.label == without.label)
        #expect(with.referenceRadiusInvAngstrom == without.referenceRadiusInvAngstrom)
    }

    /// A crystal with no forbidden shell on the detector (fcc Al [001]) is assigned as without the argument.
    @Test func anEmptyListChangesNothing() throws {
        let radii = LinearScaleAssignmentTests.positions(LinearScaleAssignmentTests.plantedRadii)
        let plain = try #require(try LinearScaleAssignmentTests.assign(radii))
        let explicit = try #require(try LinearScale.assign(radiiPerPosition: radii, shellLengths: LinearScaleAssignmentTests.al,
                                                           labels: LinearScaleAssignmentTests.alLabels, crystalName: "Al",
                                                           doubleDiffractionShells: []))
        #expect(plain == explicit)
    }

    /// The real Si-SiGe scan (`downsample_Si_SiGe_exp.h5`, mac4DSTEM#274; its ring radii as the app measured them, kept outside the
    /// repo in `~/Library/Caches/DiffractoKit/real-data/si-sige-radii.json`). With the forbidden shells the pooled reading names {111}
    /// at 15.6 px; without them it starts at {220} or {400}. The per-position reading still splits {111} (its equivalents spread about
    /// 5 % at one position, over the 3 % band), so the whole assignment is still refused: the 3 % guard is kept as it is.
    @Test func realSiSiGePooledReadingIsTheFirstShell() throws {
        let url = FileManager.default.homeDirectoryForCurrentUser.appendingPathComponent("Library/Caches/DiffractoKit/real-data/si-sige-radii.json")
        guard let data = try? Data(contentsOf: url),
              let json = try JSONSerialization.jsonObject(with: data) as? [String: Any],
              let rows = json["radii"] as? [[NSNumber]] else { return }   // not on this machine
        let radii = rows.map { $0.map(\.doubleValue) }
        let pooled = LinearScale.ringsPooled(radii, minimumRadiusPixels: 2)
        func name(_ k: Int) -> String { Self.siLabels[k] }
        let plain = try? LinearScale.assign(rings: pooled, shellLengths: Self.si, extra: [], name: name, crystalName: "Si", priorScale: nil)
        #expect(plain.flatMap { $0 }?.label != "{111}", "without the forbidden shells the pooled reading misses {111}")
        let found = try #require(try LinearScale.assign(rings: pooled, shellLengths: Self.si, extra: Self.siDoubleDiffraction,
                                                        name: name, crystalName: "Si", priorScale: nil))
        #expect(found.label == "{111}")
        #expect(abs(found.observedRadiiPixels[0] - 15.6) < 0.1)
        #expect(throws: LinearScaleAssignmentError.self) { _ = try Self.assign(radii, doubleDiffraction: Self.siDoubleDiffraction) }
    }
}
