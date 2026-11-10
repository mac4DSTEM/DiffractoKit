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

    /// Two rings at 1 : √2 fit {A}, {B} and equally {B}, {C}: the same prediction at two scales, which no radius of these rings separates.
    /// The lowest shell is assigned and the others are named (it was a refusal before mac4DSTEM#110).
    @Test func theSamePredictionAtTwoScalesTakesTheLowerShellAndNamesTheOthers() throws {
        let shells = [1.0, 2.0.squareRoot(), 2.0, 8.0.squareRoot()]
        let radii = Self.positions([40, 40 * 2.0.squareRoot()])
        let r = try #require(try Self.assign(radii, shells: shells, labels: ["{A}", "{B}", "{C}", "{D}"]))
        #expect(r.shellIndex == 0 && r.label == "{A}" && r.matchedRings == 2 && r.equallyFitting == ["{B}", "{C}"])
        let al = try #require(try Self.assign(Self.positions([Self.plantedRadii[0], Self.plantedRadii[1]])))
        #expect(al.label == "{200}" && al.equallyFitting == ["{220}"])
        #expect(try #require(try Self.assign(Self.positions(Self.plantedRadii))).equallyFitting.isEmpty)
    }

    /// Two assignments with DIFFERENT predictions that both fit are still refused and named: Al rings at 1 : 1.645 are {111}, {220}
    /// (1.633) or {200}, {311} (1.658), and a better measurement could tell them apart.
    @Test func twoEquallyGoodAssignmentsWithDifferentPredictionsAreRefusedAndNamed() {
        #expect(throws: LinearScaleAssignmentError.ambiguous(first: "{111}", second: "{200}")) {
            try Self.assign(Self.positions([40, 40 * 1.645]))
        }
        let sentence = LinearScaleAssignmentError.ambiguous(first: "{A}", second: "{B}").localizedDescription
        #expect(sentence.hasPrefix("Two assignments fit equally well ({A} and {B} first)"))
    }

    /// mac4DSTEM#110: on a binned, elliptical Al [001] pattern the four {200} spots spread 6 % in radius and the four {220} spots 7 %.
    /// A ring is the run of peaks each within 3 % of the one before, so each family stays one ring. Clustered from the ring's first
    /// peak, the outer {200} spots became a ring at ratio 1.04 and the pattern read {311}, {222}.
    @Test func aRingThatSpreadsMoreThanTheBandStaysOneRing() throws {
        let radii = [[Double]](repeating: [17.0, 17.3, 17.8, 18.0, 24.6, 25.0, 25.7, 26.3], count: 25)
        let r = try #require(try Self.assign(radii))
        #expect(r.label == "{200}" && r.observedRings == 2 && r.matchedRings == 2 && r.equallyFitting == ["{220}"])
        #expect(abs(r.observedRadiiPixels[0] - 17.525) < 1e-9 && abs(r.observedRadiiPixels[1] - 25.4) < 1e-9)
        let both = try #require(try LinearScale.estimateAssigned(radiiPerPosition: radii, shellLengths: Self.al, labels: Self.alLabels))
        #expect(abs(both.estimate.invAngstromPerPixel - Self.al[1] / 17.525) < 1e-12)
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
    /// two starting shells at the same ratio, so {200} is taken and {220} named; at 25 % the third ring is read and decides. A fourth ring is
    /// never read.
    @Test func theRingCountAndPositionShareArePinned() throws {
        let rings = [Self.plantedRadii[0], Self.plantedRadii[1]]
        func sparse(share: Double) -> [[Double]] {
            (0 ..< 100).map { i in
                let peaks = Self.positions(rings, count: 1)[0]
                return Double(i) < share * 100 ? peaks + [Self.plantedRadii[2], Self.plantedRadii[2], Self.plantedRadii[2], Self.plantedRadii[2]] : peaks
            }
        }
        let unread = try #require(try Self.assign(sparse(share: 0.15)))
        #expect(unread.observedRings == 2 && unread.shellIndex == 1 && unread.equallyFitting == ["{220}"])
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

    static let cubicLabels = ["{111}", "{200}", "{220}", "{311}", "{222}", "{400}", "{331}", "{420}"]
    static func fcc(_ a: Double) -> [Double] { [3.0, 4.0, 8.0, 11.0, 12.0, 16.0, 19.0, 20.0].map { $0.squareRoot() / a } }

    /// A polycrystal's broad first ring (mac4DSTEM#110, the py4DSTEM Au simulation): {111} spread over 19.5 to 23 px beside {200} near 25 px, then
    /// {220} and {311}. Read position by position the broad ring is two "rings" a few percent apart, which with {200} fit {311}, {222}, {400} on all
    /// three (Δq 1.9 times too large: what `assign` returned before the readings had to agree). Read over the whole scan it is one ring and the
    /// sequence is {111}. The two scales differ, so nothing is assigned. Mutation: return the per-position reading alone and {311} comes back.
    @Test func aBroadFirstRingReadTwoWaysIsRefusedNotAssignedToThreeEleven() throws {
        let au = Self.fcc(4.0782), deltaQ = 0.0195
        var state: UInt64 = 42
        func uniform() -> Double {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return Double(state >> 11) / Double(1 << 53) * 2 - 1
        }
        let r200 = au[1] / deltaQ, r220 = au[2] / deltaQ, r311 = au[3] / deltaQ
        let radii: [[Double]] = (0 ..< 400).map { i in
            var peaks = (0 ..< (3 + i % 2)).map { _ in 21.25 + 1.75 * uniform() }
            peaks += (0 ..< 2).map { _ in r200 + 0.4 * uniform() }
            peaks += [r220 + 0.5 * uniform(), r311 + 0.5 * uniform()]
            return peaks
        }
        // the two readings on their own: {311} at twice the scale, and {111}
        let perPosition = LinearScale.ringsPerPosition(radii, minimumRadiusPixels: 2), pooled = LinearScale.ringsPooled(radii, minimumRadiusPixels: 2)
        #expect(perPosition.count == 3 && perPosition[1] / perPosition[0] < 1.10, "the broad ring read as two: \(perPosition)")
        #expect(pooled.count == 3 && abs(pooled[0] - 21.25) < 0.2 && abs(pooled[1] - r200) < 0.3 && abs(pooled[2] - r220) < 0.3, "\(pooled)")
        do {
            let r = try LinearScale.assign(radiiPerPosition: radii, shellLengths: au, labels: Self.cubicLabels, crystalName: "Au")
            Issue.record("assigned \(r?.label ?? "nil"), rings \(r?.observedRadiiPixels ?? [])")
        } catch let error as LinearScaleAssignmentError {
            guard case .readingsDisagree(let a, let b) = error else { Issue.record("\(error)"); return }
            #expect(a.hasPrefix("{311} at ") && b.hasPrefix("{111} at "), "\(a), \(b)")
            #expect(error.errorDescription?.contains("read over the whole scan it is {111}") == true)
        }
    }

    /// Grains of two orientations in one scan (the Al-Mg-Si demo's shape): 70 % of the positions show Al [001] ({200}, {220}, {400}), 30 % show
    /// [011] ({111}, {200}, {220}). Position by position the first ring is {200}; over the whole scan the rings are {111}, {200}, {220}. The labels
    /// differ and the scale is the same, so the per-position assignment stands.
    @Test func twoGrainOrientationsGiveTwoLabelsAtOneScaleAndAreAssigned() throws {
        let al = Self.fcc(Self.a), deltaQ = 0.0120
        let radii: [[Double]] = (0 ..< 100).map { i in
            (i % 10 < 7 ? [al[1], al[2], al[5]] : [al[0], al[1], al[2]]).flatMap { g in [g / deltaQ, g / deltaQ] }
        }
        let pooled = LinearScale.ringsPooled(radii, minimumRadiusPixels: 2)
        #expect(abs(pooled[0] - al[0] / deltaQ) < 1e-9, "the pooled first ring is {111}: \(pooled)")
        let r = try #require(try LinearScale.assign(radiiPerPosition: radii, shellLengths: al, labels: Self.cubicLabels, crystalName: "Al"))
        #expect(r.label == "{200}")
        #expect(abs(r.referenceRadiusInvAngstrom / r.observedRadiiPixels[0] - deltaQ) <= 1e-9 * deltaQ)
    }

    /// Tight polycrystal rings of three lattices, every shell at every position: both readings see the same rings and the first shell is assigned.
    @Test func tightPolycrystalRingsOfFccBccAndDiamondAreAssignedToTheFirstShell() throws {
        let bcc = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0].map { $0.squareRoot() / 2.8665 }
        let diamond = [3.0, 8.0, 11.0, 16.0, 19.0, 24.0, 27.0, 32.0].map { $0.squareRoot() / 5.4309 }
        for (shells, deltaQ) in [(Self.fcc(Self.a), 0.027), (bcc, 0.035), (diamond, 0.020)] {
            let radii = Self.positions(shells.prefix(4).map { $0 / deltaQ }, noise: 0.1, seed: 3)
            let r = try #require(try LinearScale.assign(radiiPerPosition: radii, shellLengths: shells, crystalName: "x"))
            #expect(r.shellIndex == 0 && r.matchedRings == 3)
            #expect(abs(r.referenceRadiusInvAngstrom / r.observedRadiiPixels[0] / deltaQ - 1) < 0.005)
        }
    }
}
