//
//  IndexingDisplayTests.swift
//  `IndexingDisplay` (ported from macTEM's IndexingDisplay.swift; the app source is the oracle, no Python reference).
//  Displayed numbers are pinned (Fixtures/indexing-display-pins.json) by tools/mactem-pieces/indexing-display-pins.py, a
//  pure-Python transcription of the app source's formulas, NOT this Swift. The decision itself is the kit's RingIndexing (P9);
//  these tests check the derived sigma_q, delta, z, margin, odds and the text, the refusal names, and that the decision is untouched.
//

import Foundation
import Testing
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

private struct DisplayPins: Decodable {
    struct Case: Decodable {
        let q, tolerancePercent, sigmaQ, g, z, deltaPercent, relativeDeltaD, dMeasured, dListed, nextG, margin: Double
        let odds: Double?
        let oddsInfinite: Bool
        let oddsText: String
    }
    let families: [Double]
    let cases: [Case]
}

struct IndexingDisplayTests {

    private static func pins() throws -> DisplayPins {
        let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/indexing-display-pins.json")
        return try JSONDecoder().decode(DisplayPins.self, from: Data(contentsOf: url))
    }

    private func close(_ a: Double?, _ b: Double, _ rel: Double = 1e-12) -> Bool {
        guard let a else { return false }
        return abs(a - b) <= rel * max(1, abs(b))
    }

    @Test func displayedNumbersMatchThePythonTranscription() throws {
        let p = try Self.pins()
        let list = try IndexPresets.dyadic(p.families)
        for c in p.cases {
            let d = IndexingDisplay.index(q: c.q, tolerancePercent: c.tolerancePercent, families: list, kMax: 4)
            #expect(d.verdict == "accepted", "q \(c.q) tol \(c.tolerancePercent): \(d.verdict)")
            #expect(d.isAccepted && d.refusal == nil && d.refusalText == nil)
            #expect(close(d.sigmaQ, c.sigmaQ) && d.q == c.q)
            #expect(d.family?.g == c.g)
            #expect(close(d.z, c.z) && close(d.deltaPercent, c.deltaPercent) && close(d.relativeDeltaD, c.relativeDeltaD))
            #expect(close(d.dMeasured, c.dMeasured) && close(d.dListed, c.dListed))
            #expect(close(d.margin, c.margin, 1e-11))
            if c.oddsInfinite { #expect(d.odds.map { $0 > 1e6 } == true) } else { #expect(close(d.odds, c.odds!, 1e-11)) }
            #expect(d.oddsText == c.oddsText)
            #expect(d.oddsBeyondCap == (c.oddsText == DisplayedIndexing.oddsCapLabel))
            #expect(d.validation == RingIndexing.validation)
        }
    }

    @Test func theDecisionIsTheKitsUntouched() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2, 3])
        let d = IndexingDisplay.index(q: 1.1, tolerancePercent: 10, families: list, kMax: 4)
        let direct = try RingIndexing.index(q: 1.1, sigma: 1.1 * 10 / 100, in: list, kMax: 4)
        #expect(d.assignment == direct)
        #expect(DisplayedIndexing.oddsCap == 1e6 && DisplayedIndexing.oddsCapLabel == "> 10⁶ : 1")
    }

    @Test func refusalsKeepTheKitsNamesAndNoNumbers() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2, 3])
        func verdict(_ q: Double, _ tol: Double, _ k: Double = 3) -> DisplayedIndexing {
            IndexingDisplay.index(q: q, tolerancePercent: tol, families: list, kMax: 4, coverage: k)
        }
        let ambiguous = verdict(1.25, 10, 2)                         // sigma 0.125: z = 2 and 2 at k = 2
        #expect(ambiguous.verdict == "ambiguous" && !ambiguous.isAccepted)
        #expect(ambiguous.family == nil && ambiguous.deltaPercent == nil && ambiguous.z == nil && ambiguous.margin == nil)
        #expect(ambiguous.odds == nil && ambiguous.oddsText == "—" && !ambiguous.oddsBeyondCap)
        #expect(ambiguous.refusalText?.hasPrefix("More than one family fits") == true)
        #expect(verdict(1.3, 1).verdict == "noneWithinTolerance")
        #expect(verdict(3.9, 10).verdict == "beyondKMax")
        for bad in [0.0, -1, .nan] {
            let d = verdict(bad, 10)
            #expect(d.verdict == "invalidInput" && d.refusalText?.hasPrefix("Cannot index:") == true)
        }
        let zeroTol = verdict(1.1, 0)
        #expect(zeroTol.verdict == "invalidInput" && zeroTol.sigmaQ == 0)
        let nanTol = verdict(1.1, .nan)
        #expect(nanTol.verdict == "invalidInput" && nanTol.sigmaQ.isNaN)
    }

    @Test func theKinematicListOverloadUsesTheListsKMax() throws {
        let list = try IndexPresets.gold.kinematicReflections(kMax: 1.2)
        let g111 = 3.0.squareRoot() / 4.0782
        let d = IndexingDisplay.index(q: g111 * 1.001, tolerancePercent: 0.5, in: list)
        #expect(d.verdict == "accepted")
        #expect(d.family?.name == "{1 1 1}")
        #expect(close(d.deltaPercent, 100 * (1 / 1.001 - 1), 1e-9))
        #expect(IndexingDisplay.index(q: 1.19, tolerancePercent: 5, in: list).verdict != "accepted")   // reach passes k_max = 1.2
    }
}
