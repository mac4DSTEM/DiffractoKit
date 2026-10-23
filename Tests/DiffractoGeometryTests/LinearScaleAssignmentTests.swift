import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoGeometry

/// `LinearScale.assign` against planted truth (calibration design 2026-10-22, section 6). No Python reference exists; every number is a closed form.
@Suite("LinearScale.assign")
struct LinearScaleAssignmentTests {
    static let a = 4.0495
    /// Al's distinct shells, 1/d in Å⁻¹: 111, 200, 220, 311, 222, 400.
    static let al: [Double] = [3.0, 4.0, 8.0, 11.0, 12.0, 16.0].map { $0.squareRoot() / a }
    static let alLabels = ["{111}", "{200}", "{220}", "{311}", "{222}", "{400}"]
    static let deltaQ = 0.0100
    /// Radii (px) of the Al [001] rings {200}, {220}, {400} at Δq = 0.0100.
    static var plantedRadii: [Double] { [al[1], al[2], al[5]].map { $0 / deltaQ } }

    static func positions(_ rings: [Double], count: Int = 25, noise: Double = 0, seed: UInt64 = 1, descan: Double = 0) -> [[Double]] {
        var state = seed
        func uniform() -> Double {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return Double(state >> 11) / Double(1 << 53) * 2 - 1
        }
        return (0 ..< count).map { _ in
            let shift = descan * uniform()
            return rings.flatMap { r in (0 ..< 4).map { k in r + shift * (k % 2 == 0 ? 1 : -1) + noise * uniform() } }
        }
    }

    static func assign(_ radii: [[Double]], shells: [Double] = al, labels: [String]? = alLabels, prior: Double? = nil) throws -> LinearScaleAssignment? {
        try LinearScale.assign(radiiPerPosition: radii, shellLengths: shells, labels: labels, crystalName: "Al", priorScale: prior)
    }

    @Test func al001FirstRingIsTwoHundredAndDeltaQRecovers() throws {
        let r = try #require(try Self.assign(Self.positions(Self.plantedRadii)))
        #expect(r.shellIndex == 1 && r.label == "{200}" && r.matchedRings == 3 && r.observedRings == 3 && r.isChecked && r.validation == "none")
        #expect(abs(r.referenceRadiusInvAngstrom / r.observedRadiiPixels[0] - Self.deltaQ) <= 1e-9 * Self.deltaQ)
        #expect(abs(try #require(r.nextShellInvAngstrom) - Self.al[2]) < 1e-12)
        // today's assumption (shell 0 as the first ring) reads Δq 13.4 % low: the defect this fixes
        let assumed = try #require(LinearScale.estimate(radiiPerPosition: Self.positions(Self.plantedRadii), referenceRadiusInvAngstrom: Self.al[0],
                                                        secondShellRadiusInvAngstrom: Self.al[1]))
        #expect(abs(assumed.invAngstromPerPixel / Self.deltaQ - 1 + 0.134) < 0.002)
        let both = try #require(try LinearScale.estimateAssigned(radiiPerPosition: Self.positions(Self.plantedRadii), shellLengths: Self.al, labels: Self.alLabels))
        #expect(abs(both.estimate.invAngstromPerPixel - Self.deltaQ) <= 1e-9 * Self.deltaQ && both.assignment.shellIndex == 1)
    }

    @Test func al001WithNoiseAndDescan() throws {
        for radii in [Self.positions(Self.plantedRadii, noise: 0.3, seed: 7), Self.positions(Self.plantedRadii, seed: 11, descan: 1.0)] {
            let r = try #require(try Self.assign(radii))
            #expect(r.shellIndex == 1 && r.matchedRings == 3)
            #expect(abs(r.referenceRadiusInvAngstrom / r.observedRadiiPixels[0] - Self.deltaQ) <= 0.005 * Self.deltaQ)
        }
    }

    @Test func selfSimilarRingsAreNotAmbiguousBecauseHigherShellsAreOutsideTheWindow() throws {
        let longList = [3.0, 4, 8, 11, 12, 16, 19, 20, 24, 27, 32, 35, 36, 40].map { $0.squareRoot() / Self.a }
        let r = try #require(try Self.assign(Self.positions(Self.plantedRadii), shells: longList, labels: nil))
        #expect(r.shellIndex == 1 && r.matchedRings == 3)
    }

    @Test func ringsFromTheFirstShellStayAssignedToIt() throws {
        let rings = [Self.al[0], Self.al[1], Self.al[2]].map { $0 / 0.0125 }
        let r = try #require(try Self.assign(Self.positions(rings)))
        #expect(r.shellIndex == 0 && r.matchedRings == 3)
    }

    @Test func twoEquallyGoodAssignmentsAreRefusedAndNamed() {
        let shells = [1.0, 2.0.squareRoot(), 2.0, 8.0.squareRoot()]
        let radii = Self.positions([40, 40 * 2.0.squareRoot()])
        #expect(throws: LinearScaleAssignmentError.ambiguous(first: "{A}", second: "{B}")) {
            try Self.assign(radii, shells: shells, labels: ["{A}", "{B}", "{C}", "{D}"])
        }
        let sentence = LinearScaleAssignmentError.ambiguous(first: "{A}", second: "{B}").localizedDescription
        #expect(sentence.hasPrefix("Two assignments fit equally well ({A} and {B} first)"))
    }

    @Test func ringsThatMatchNoStartingShellAreRefusedWithTheRatios() throws {
        do { _ = try Self.assign(Self.positions([50, 65])); Issue.record("expected a refusal") } catch let e as LinearScaleAssignmentError {
            let s = e.localizedDescription
            #expect(s.hasPrefix("The rings do not match any starting shell of Al: the first ring gives ratios 1.30"), Comment(rawValue: s))
            #expect(s.contains("({111} first)") && s.contains("({200} first)") && s.contains("1.63"), Comment(rawValue: s))
        }
    }

    @Test func aSingleRingIsAssumedToBeTheFirstShellAndNotChecked() throws {
        let r = try #require(try Self.assign(Self.positions([Self.al[0] / Self.deltaQ])))
        #expect(r.shellIndex == 0 && !r.isChecked && r.observedRings == 1)
    }

    @Test func aPriorScaleBreaksOnlyTheSingleRingTie() throws {
        let one = Self.positions([Self.al[1] / Self.deltaQ])
        #expect(try #require(try Self.assign(one, prior: 0.0100)).shellIndex == 1)
        #expect(try #require(try Self.assign(one, prior: 0.0090)).shellIndex == 0)
        #expect(try #require(try Self.assign(one, prior: 5)).shellIndex == 0)
        #expect(try #require(try Self.assign(Self.positions(Self.plantedRadii), prior: 0.00866)).shellIndex == 1)   // the sequence wins
    }

    @Test func noUsableRadiusGivesNil() throws {
        #expect(try Self.assign([[], [1, .nan]]) == nil)
        #expect(try LinearScale.assign(radiiPerPosition: [[50]], shellLengths: []) == nil)
        #expect(try LinearScale.estimateAssigned(radiiPerPosition: [[]], shellLengths: Self.al) == nil)
    }

    /// The sweep that pins K = 3 rings and the 20 % share. {400} seen at 15 % of the positions is not read: the two rings left ({200}, {220}) fit
    /// two starting shells and the call refuses as ambiguous; at 25 % the third ring is read and decides. A fourth ring is never read.
    @Test func theRingCountAndPositionShareArePinned() throws {
        let rings = [Self.plantedRadii[0], Self.plantedRadii[1]]
        func sparse(share: Double) -> [[Double]] {
            (0 ..< 100).map { i in
                let peaks = Self.positions(rings, count: 1)[0]
                return Double(i) < share * 100 ? peaks + [Self.plantedRadii[2], Self.plantedRadii[2], Self.plantedRadii[2], Self.plantedRadii[2]] : peaks
            }
        }
        #expect(throws: LinearScaleAssignmentError.self) { try Self.assign(sparse(share: 0.15)) }
        let kept = try #require(try Self.assign(sparse(share: 0.25)))
        #expect(kept.observedRings == 3 && kept.shellIndex == 1)
        // a fourth ring (here the Al {422}-like radius) at every position: still three rings read
        let four = Self.positions(Self.plantedRadii + [Self.al[5] / Self.deltaQ * 1.22])
        #expect(try #require(try Self.assign(four)).observedRings == 3)
    }

    @Test func labelsAreWrittenCompactWhenEveryIndexIsOneDigit() {
        #expect(LinearScale.compactShellLabel("{2 0 0}") == "{200}")
        #expect(LinearScale.compactShellLabel("{1 0 -1 1}") == "{1 0 -1 1}")
        #expect(LinearScale.compactShellLabel("{3 3 3} + {5 1 1}") == "{3 3 3} + {5 1 1}")
        #expect(LinearScale.compactShellLabel("{10 0 0}") == "{10 0 0}")
    }
}
