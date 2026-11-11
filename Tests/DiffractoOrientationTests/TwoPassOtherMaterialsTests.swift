//
//  TwoPassOtherMaterialsTests.swift
//  The two passes (DiffractoKit#81) on matrix and precipitate pairs other than aluminium's: nothing in them may be tuned to one
//  material. Planted truth, the same scan for every pair:
//    - 24 positions with the full pattern of the precipitate's richest variants over the matrix pattern,
//    - 12 positions with three of the precipitate's own spots over the matrix pattern,
//    -  6 positions with a rich pattern turned 45 degrees about the beam (a crystal that does not keep the relationship),
//    -  4 positions with ONE precipitate spot over the matrix pattern (an edge: the firm rule calls it matrix),
//    - 12 positions of the matrix alone.
//  The matrix removal is emulated exactly: the residual is every precipitate spot farther than the removal's tolerance from
//  each matrix spot, so a precipitate that shares spots with its matrix loses them, as it does in a measured scan.
//  Pairs: hexagonal matrix with a cubic precipitate (alpha-Ti, beta-Ti, Burgers); cubic matrix with an ordered precipitate whose
//  fundamental spots coincide with the matrix's (Ni, gamma' Ni3Al, cube on cube); bcc matrix with an orthorhombic precipitate
//  (alpha-Fe, cementite Fe3C, Bagaryatsky); and the monoclinic precipitate in a cubic matrix the other suites use (Al, beta'').
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct TwoPassOtherMaterialsTests {
    static let q = 0.012
    static let origin: Float = 63.5
    static let tolerance = 0.02

    struct Pair {
        let name: String
        let matrix: CrystalStructure
        let precipitate: CrystalStructure
        /// Two parallel pairs that fix the relationship: precipitate element, matrix element.
        let first: OrientationRelationship
        let second: OrientationRelationship
        let reach: Double
    }

    struct Outcome: CustomStringConvertible {
        var name = ""
        var plantedVariants = 0, learnedRelationships = 0, learnedVariants = 0
        var relationErrorDeg = Double.nan
        var rich = 0, richRight = 0, weak = 0, weakRight = 0, weakFreeSearch = 0
        var turned = 0, turnedFlagged = 0, single = 0, singleCalled = 0, matrix = 0, matrixCalled = 0
        var richSpots = 0, ownSpotsOfRich = 0
        var besideCalled = false, besideSupported = false
        var matrixIsMatrix = 0, singleIsMatrix = 0
        var description: String {
            String(format: "%@: planted variants %d | learned %d relationship(s), %d variants, %.2f deg off | rich %d/%d (pattern %d spots, %d its own) | weak %d/%d (free search alone %d) | turned flagged %d/%d | one-spot called %d/%d | matrix called %d/%d",
                   name, plantedVariants, learnedRelationships, learnedVariants, relationErrorDeg, richRight, rich, richSpots, ownSpotsOfRich,
                   weakRight, weak, weakFreeSearch, turnedFlagged, turned, singleCalled, single, matrixCalled, matrix)
        }
    }

    static func cementite() -> CrystalStructure {
        // Fe3C, Pnma, a 5.0896, b 6.7443, c 4.5248 Å (Fasiska and Jeffrey 1965): Fe 8d, Fe 4c, C 4c.
        func general(_ p: SIMD3<Double>) -> [SIMD3<Double>] {
            let (x, y, z) = (p.x, p.y, p.z)
            let half: [SIMD3<Double>] = [[x, y, z], [-x + 0.5, -y, z + 0.5], [-x, y + 0.5, -z], [x + 0.5, -y + 0.5, -z + 0.5]]
            return (half + half.map { -$0 }).map { v in SIMD3(v.x - floor(v.x), v.y - floor(v.y), v.z - floor(v.z)) }
        }
        func distinct(_ list: [SIMD3<Double>]) -> [SIMD3<Double>] {
            var out: [SIMD3<Double>] = []
            for v in list where !out.contains(where: { o in
                let d = v - o
                return simd_length(SIMD3(d.x - d.x.rounded(), d.y - d.y.rounded(), d.z - d.z.rounded())) < 1e-4
            }) { out.append(v) }
            return out
        }
        let fe = distinct(general([0.1816, 0.0666, 0.3374])) + distinct(general([0.0367, 0.25, 0.8402]))
        let c = distinct(general([0.8764, 0.25, 0.4426]))
        return CrystalPresets.structure(a: 5.0896, b: 6.7443, c: 4.5248,
                                        sites: fe.map { AtomSite(z: 26, fractional: $0) } + c.map { AtomSite(z: 6, fractional: $0) })
    }

    static var pairs: [Pair] {
        [
            // Burgers: (0001)alpha ∥ (110)beta, [2-1-10]alpha ∥ [1-11]beta.
            Pair(name: "alpha-Ti (hcp) + beta-Ti (bcc), Burgers",
                 matrix: CrystalPresets.hcp(a: 2.951, c: 4.684, z: 22), precipitate: CrystalPresets.bcc(a: 3.306, z: 22),
                 first: OrientationRelationship(candidate: .plane(SIMD3(1, 1, 0)), matrix: .plane(SIMD3(0, 0, 1))),
                 second: OrientationRelationship(candidate: .direction(SIMD3(1, -1, 1)), matrix: .direction(SIMD3(1, 0, 0))),
                 reach: 1.2),
            // Cube on cube; gamma' is L1_2 (Al at the corner, Ni on the faces): its fundamental spots are the matrix's.
            Pair(name: "Ni (fcc) + gamma' Ni3Al (L1_2), cube on cube",
                 matrix: CrystalPresets.fcc(a: 3.524, z: 28),
                 precipitate: CrystalPresets.structure(a: 3.572, b: 3.572, c: 3.572, sites: [
                    AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 28, fractional: [0.5, 0.5, 0]),
                    AtomSite(z: 28, fractional: [0.5, 0, 0.5]), AtomSite(z: 28, fractional: [0, 0.5, 0.5])]),
                 first: OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(1, 0, 0))),
                 second: OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 1, 0))),
                 reach: 0.75),
            // Bagaryatsky: [100]c ∥ [1-10]alpha, [010]c ∥ [111]alpha.
            Pair(name: "alpha-Fe (bcc) + cementite Fe3C (orthorhombic), Bagaryatsky",
                 matrix: CrystalPresets.bcc(a: 2.866, z: 26), precipitate: cementite(),
                 first: OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(1, -1, 0))),
                 second: OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(1, 1, 1))),
                 reach: 0.75),
            Pair(name: "Al (fcc) + beta'' (monoclinic), Andersen",
                 matrix: CrystalPresets.aluminum, precipitate: CrystalPresets.betaDoublePrime,
                 first: OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
                 second: OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(3, 1, 0))),
                 reach: 0.75),
        ]
    }

    static func run(_ pair: Pair, usesNeighbourSupport: Bool = false) throws -> Outcome {
        var out = Outcome(name: pair.name)
        let relation = try OrientationRelationship.rotation(pair.first, pair.second, precipitate: pair.precipitate, matrix: pair.matrix)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(pair.matrix.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(pair.precipitate.laueClass))
        let planted = LearnedOrientationRelationship(relation: relation, positionCount: 0, spreadRad: 0, matrixSymmetry: ms, precipitateSymmetry: ps)
        out.plantedVariants = planted.variants.count
        let candidate = try PhaseOrientationCandidate(name: "precipitate", crystal: pair.precipitate, kMax: pair.reach)
        let matrixCandidate = try PhaseOrientationCandidate(name: "matrix", crystal: pair.matrix, kMax: pair.reach)
        // The matrix along its c axis (cubic [001], hexagonal [0001]), turned 12 degrees in the plane.
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(12 * Double.pi / 180))
        func inReach(_ spots: [SIMD2<Double>]) -> [SIMD2<Double>] {
            spots.filter { simd_length($0) > 0.15 && simd_length($0) < pair.reach }.sorted { simd_length($0) < simd_length($1) }
        }
        let matrixSpots = inReach(PhaseVectorMatcher.matrixSpots(reflections: matrixCandidate.reflections, orientation: m, excitationSlabInvAngstrom: 0.02))
        func spots(_ variant: simd_double3x3, turnedDeg: Double = 0) -> [SIMD2<Double>] {
            let c = cos(turnedDeg * .pi / 180), s = sin(turnedDeg * .pi / 180)
            return inReach(PhaseVectorMatcher.matrixSpots(reflections: candidate.reflections, orientation: variant.transpose * m, excitationSlabInvAngstrom: 0.02))
                .map { SIMD2(c * $0.x - s * $0.y, s * $0.x + c * $0.y) }
        }
        func own(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] {
            list.filter { p in !matrixSpots.contains { simd_length($0 - p) <= tolerance } }
        }
        func pixels(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] { list.map { SIMD2(Double(origin) + $0.x / q, Double(origin) + $0.y / q) } }
        // The richest variants: those with the most spots of their own.
        let byRichness = planted.variants.sorted { own(spots($0)).count > own(spots($1)).count }
        let most = own(spots(byRichness[0])).count
        let rich = byRichness.filter { own(spots($0)).count == most }
        let others = byRichness.filter { v in own(spots(v)).count >= 3 && !rich.contains { $0 == v } }
        out.richSpots = spots(rich[0]).count
        out.ownSpotsOfRich = most

        var scan: [[SIMD2<Double>]] = [], residual: [[SIMD2<Double>]] = [], kind: [String] = []
        func add(_ precipitateSpots: [SIMD2<Double>], _ k: String) {
            scan.append(pixels(matrixSpots + precipitateSpots))
            residual.append(pixels(own(precipitateSpots)))
            kind.append(k)
        }
        for k in 0..<24 { add(spots(rich[k % rich.count]), "rich") }
        for k in 0..<12 {
            let source = others.isEmpty ? rich : others
            add(Array(own(spots(source[k % source.count])).prefix(3)), "weak")
        }
        for _ in 0..<6 { add(spots(rich[0], turnedDeg: 45), "turned") }
        for k in 0..<4 { add(Array(own(spots(rich[k % rich.count])).prefix(1)), "single") }
        for _ in 0..<12 { add([], "matrix") }
        // The tail, for the neighbour rule (the scan is one row): a rich position, one spot of the same variant beside it, two
        // matrix positions, the same one spot alone, a matrix position.
        add(spots(rich[0]), "rich")
        add(Array(own(spots(rich[0])).prefix(1)), "single beside")
        add([], "matrix"); add([], "matrix")
        add(Array(own(spots(rich[0])).prefix(1)), "single")
        add([], "matrix")

        let result = try #require(PhaseVectorMatcher.classifyInTwoPasses(
            residual: MatrixResidualScan(scan: residual, removedCounts: [Int](repeating: matrixSpots.count, count: scan.count),
                                         inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count)),
            matrixOrientations: [simd_double3x3?](repeating: m, count: scan.count),
            phases: [TwoPassPhase(candidate: candidate, symmetry: ps)], matrixSymmetry: ms,
            settings: PhaseOrientationSearchSettings(reachInvAngstrom: pair.reach),
            originX: origin, originY: origin, invAngstromPerPixel: q,
            matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scan.count,
            usesNeighbourSupport: usesNeighbourSupport, minimumPositions: 10))
        out.learnedRelationships = result.relationships[0].count
        if let best = result.relationships[0].max(by: { $0.positionCount < $1.positionCount }) {
            out.learnedVariants = best.variants.count
            out.relationErrorDeg = PhaseVectorMatcher.relationDistance(best.relation, relation, matrixSymmetry: ms, precipitateSymmetry: ps) * 180 / .pi
        }
        for (i, k) in kind.enumerated() {
            let v = result.verdicts[i], called = v.verdict == .indexed
            switch k {
            case "rich": out.rich += 1; if called { out.richRight += 1 }
            case "weak":
                out.weak += 1
                if called { out.weakRight += 1 }
                if result.freeSearch[i].verdict == .indexed { out.weakFreeSearch += 1 }
            case "turned":
                out.turned += 1
                // Flagged: not called at the relationship, or called by the free search far from it.
                let angle = result.relationAnglesRad[i] * 180 / .pi
                if !called || angle > 10 { out.turnedFlagged += 1 }
            case "single": out.single += 1; if called { out.singleCalled += 1 }; if v.verdict == .matrix { out.singleIsMatrix += 1 }
            case "single beside": out.besideCalled = called; out.besideSupported = result.neighbourSupported[i]
            default: out.matrix += 1; if called { out.matrixCalled += 1 }; if v.verdict == .matrix { out.matrixIsMatrix += 1 }
            }
        }
        return out
    }

    /// What each pair is held to, as measured 2026-11-11.
    /// - alpha-Ti + beta-Ti: every position right. The relationship comes back 2.65 degrees off with twelve variants for six:
    ///   Burgers' variants come in pairs 5.26 degrees apart, the free search does not tell the two of a pair apart, and the
    ///   learned relation sits between them.
    /// - Ni + gamma': KNOWN GAP. Nothing is learned and no gamma' is found. The residual holds only the sixteen superlattice
    ///   spots (the fundamental ones are the matrix's and go with the removal), and the free search never proposes the true
    ///   orientation: its best has 10 of 16 spots at a score of 0.22. Asking for up to 200 proposals or flattening the
    ///   template weights (intensity power 0.1, 0.01) does not change that. The expectation pins today's behaviour, so the
    ///   test fails, and is to be turned round, the day the free search finds it.
    /// - alpha-Fe + cementite: every position right, the relationship 1.44 degrees off.
    /// - Al + beta'': every position right, the relationship exact.
    static let expected: [(learns: Bool, maximumErrorDeg: Double, variants: Int)] = [
        (true, 3, 12), (false, .nan, 0), (true, 2, 12), (true, 0.5, 12),
    ]

    /// Every pair: the relationship comes back, the rich and the weak positions are the phase, the matrix and the one-spot
    /// positions are not, and the turned crystal is flagged.
    @Test(arguments: 0..<4) func thePairIsClassifiedFromItsStructuresAlone(_ index: Int) throws {
        let o = try Self.run(Self.pairs[index]), want = Self.expected[index]
        print("OTHER MATERIALS | \(o)")
        #expect(o.matrixCalled == 0 && o.matrixIsMatrix == o.matrix, "every matrix position's verdict is .matrix, also where nothing was learned. \(o)")
        #expect(o.singleCalled == 0 && o.singleIsMatrix == o.single && !o.besideCalled, "\(o)")
        guard want.learns else {
            #expect(o.learnedRelationships == 0 && o.richRight == 0 && o.weakRight == 0, "the known gap is closed: turn this expectation round. \(o)")
            return
        }
        #expect(o.learnedRelationships >= 1, "\(o)")
        #expect(o.relationErrorDeg < want.maximumErrorDeg, "\(o)")
        #expect(o.learnedVariants == want.variants, "\(o)")
        #expect(o.richRight == o.rich, "\(o)")
        #expect(o.weakRight == o.weak, "\(o)")
        #expect(o.turnedFlagged == o.turned, "\(o)")
    }

    /// The neighbour rule: one spot of a precipitate beside a firm position of the same variant is the phase and is marked as
    /// supported; the same spot alone, and every matrix position, stay matrix.
    @Test(arguments: [0, 2, 3]) func oneSpotBesideAFirmPositionIsThePhaseWhenAskedFor(_ index: Int) throws {
        let o = try Self.run(Self.pairs[index], usesNeighbourSupport: true)
        print("OTHER MATERIALS, neighbour rule | \(o) | beside: called \(o.besideCalled), supported \(o.besideSupported)")
        #expect(o.besideCalled && o.besideSupported, "\(o)")
        #expect(o.singleCalled == 0 && o.singleIsMatrix == o.single, "\(o)")
        #expect(o.matrixCalled == 0 && o.matrixIsMatrix == o.matrix, "\(o)")
        #expect(o.richRight == o.rich && o.weakRight == o.weak, "\(o)")
    }

    /// Two relationships of one phase within the sum of their spreads become one: the mean weighted by their positions (30
    /// and 10 positions 2 degrees apart: 0.5 degrees from the larger), the positions added, the spread covering both. Two
    /// relationships farther apart than their spreads stay two.
    @Test func relationshipsWithinTheirSpreadsAreMerged() throws {
        let pair = Self.pairs[3]
        let relation = try OrientationRelationship.rotation(pair.first, pair.second, precipitate: pair.precipitate, matrix: pair.matrix)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(pair.matrix.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(pair.precipitate.laueClass))
        func turned(_ deg: Double) -> simd_double3x3 {
            simd_double3x3(simd_quatd(angle: deg * .pi / 180, axis: simd_normalize(SIMD3(1, 2, 3)))) * relation
        }
        func learned(_ r: simd_double3x3, _ count: Int, spreadDeg: Double) -> LearnedOrientationRelationship {
            LearnedOrientationRelationship(relation: r, positionCount: count, spreadRad: spreadDeg * .pi / 180, matrixSymmetry: ms, precipitateSymmetry: ps)
        }
        func distanceDeg(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
            PhaseVectorMatcher.relationDistance(a, b, matrixSymmetry: ms, precipitateSymmetry: ps) * 180 / .pi
        }
        let one = PhaseVectorMatcher.mergingCloseRelationships([learned(turned(2), 10, spreadDeg: 1.5), learned(relation, 30, spreadDeg: 1.5)],
                                                               matrixSymmetry: ms, precipitateSymmetry: ps)
        #expect(one.count == 1)
        let merged = try #require(one.first)
        #expect(merged.positionCount == 40)
        #expect(abs(distanceDeg(merged.relation, relation) - 0.5) < 0.05, "\(distanceDeg(merged.relation, relation)) deg from the larger")
        #expect(abs(distanceDeg(merged.relation, turned(2)) - 1.5) < 0.05)
        #expect(merged.spreadRad * 180 / .pi >= 2 - 1e-6)
        let two = PhaseVectorMatcher.mergingCloseRelationships([learned(relation, 30, spreadDeg: 1.5), learned(turned(4), 10, spreadDeg: 1.5)],
                                                               matrixSymmetry: ms, precipitateSymmetry: ps)
        #expect(two.count == 2)
    }
}
