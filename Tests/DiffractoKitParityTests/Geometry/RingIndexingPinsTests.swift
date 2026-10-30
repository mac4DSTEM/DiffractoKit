//
//  RingIndexingPinsTests.swift
//  The indexing decision against the pins written BEFORE its code:
//    numpy-ring-indexing-pins.json  tools/pins/ring-indexing/pins.py: NumPy from the written definition on the kit's own ring lists
//                             (DiffractoScattering's families after the cut), T1 (every family x u x s) and T4 (refusal plants).
//  Verdict equal, payload family names equal and in order, delta 1e-14, z 1e-10 max(1, z), margin 1e-9 max(1, margin),
//  reach 1e-9 max(1, reach). A failure is a failure of the port, never a tolerance to widen.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

// MARK: - The pins, decoded

private struct PinList: Decodable {
    let name: String
    let source: String
    let kmax: Double
    let names: [String]
    let g: [Double]
}

private struct PinPlant: Decodable {
    let c: Int?             // T1: list index
    let t: Int?             // T1: truth family
    let u: Double?
    let s: Double?
    let l: String?          // T4: list name
    let label: String?
    let q: Double
    let sg: Double
    let k: Double?
    let v: String           // verdict
    let b: Int?             // best
    let z: Double?
    let d: Double?          // delta
    let n: Int?             // next
    let m: Double?          // margin
    let p: [Int]?           // payload (ambiguous)
    let nr: Int?            // nearest (none)
    let nz: Double?
    let nd: Double?
    let reach: Double?
    let expect: String?
}

private struct PinHazards: Decodable {
    let minAbsZMinusK: Double
    let minAbsMarginMinusLn19: Double
    let minAbsReachMinusKmax: Double
}

private struct PinDefinitions: Decodable {
    let k: Double
    let ln19: Double
    let kMax: Double
    let u: [Double]
    let s: [Double]
}

private struct IndexPins: Decodable {
    let definitions: PinDefinitions
    let lists: [PinList]
    let t1Lists: [String]
    let t1: [PinPlant]
    let t4: [PinPlant]
    let hazards: PinHazards
}

private enum IndexPinFixtures {
    static let directory: URL = FixtureDirs.target("DiffractoGeometryTests").appendingPathComponent("Fixtures")

    static let result: Result<IndexPins, any Error> = Result {
        let data = try Data(contentsOf: directory.appendingPathComponent("numpy-ring-indexing-pins.json"))
        let d = JSONDecoder()
        d.keyDecodingStrategy = .convertFromSnakeCase
        d.nonConformingFloatDecodingStrategy = .convertFromString(positiveInfinity: "Infinity", negativeInfinity: "-Infinity",
                                                                  nan: "NaN")
        return try d.decode(IndexPins.self, from: data)
    }

    static func pins() throws -> IndexPins { try result.get() }

    /// The structure each pinned list was derived from (the kit's own builders, copied into this target).
    static func structure(_ source: String) throws -> CrystalStructure {
        switch source {
        case "Au": return IndexPresets.gold
        case "Al": return IndexPresets.aluminium
        case "Si": return IndexPresets.silicon
        case "Mg": return IndexPresets.magnesium
        case "ideal_hcp": return IndexPresets.idealHCP
        case "Weak-cut": return IndexPresets.weakCut
        case "beta_double_prime": return IndexPresets.betaDoublePrime
        default:
            Issue.record("no builder for pinned structure \(source)")
            throw IndexPinError.unknownStructure(source)
        }
    }

    /// The Swift rings of a pinned list, at its k_max.
    static func rings(_ list: PinList) throws -> [KinematicFamily] {
        try structure(list.source).kinematicReflections(kMax: list.kmax).families
    }
}

private enum IndexPinError: Error { case unknownStructure(String) }

private func within(_ x: Double, _ want: Double, _ tol: Double) -> Bool {
    abs(x - want) <= tol
}

struct RingIndexingPinsTests {

    private typealias Index = RingIndexing

    /// The pins carry their definitions (k = 3, ln 19, k_max 1.2, the five u and nine s), and the plant counts.
    @Test func testPinsCarryTheirDefinitions() throws {
        let pins = try IndexPinFixtures.pins()
        #expect(pins.definitions.k == Index.defaultCoverage)
        #expect(abs(pins.definitions.ln19 - Index.minimumMargin) <= 1e-15)
        #expect(pins.definitions.kMax == 1.2)
        #expect(pins.definitions.u == [1e-4, 1e-3, 3e-3, 1e-2, 3e-2])
        #expect(pins.definitions.s == [-2.9, -2, -1, -0.5, 0, 0.5, 1, 2, 2.9])
        #expect(pins.t1Lists == ["Au", "Al", "Si", "Mg", "Mg-ideal", "Weak-cut", "beta''"])
        #expect(pins.t1.count == (8 + 8 + 10 + 24 + 23 + 10 + 417) * 45)
        #expect(pins.t4.count >= 35)
        // the plants sit away from the three edges the two implementations could round across
        #expect(pins.hazards.minAbsZMinusK > 1e-6)
        #expect(pins.hazards.minAbsMarginMinusLn19 > 1e-6)
        #expect(pins.hazards.minAbsReachMinusKmax > 1e-9)
    }

    /// The pinned ring lists are the kit's own: the same count, the same names in the same order, g to 1e-12 relative.
    @Test func testPinnedRingListsAreTheKitsRings() throws {
        let pins = try IndexPinFixtures.pins()
        #expect(pins.lists.map(\.name) == ["Au", "Al", "Si", "Mg", "Mg-ideal", "Weak-cut", "beta''", "Au@0.45", "Au@0.3"])
        for list in pins.lists {
            let rings = try IndexPinFixtures.rings(list)
            #expect(rings.map(\.family.name) == list.names, "\(list.name)")
            #expect(rings.count == list.g.count, "\(list.name)")
            for (r, g) in zip(rings, list.g) where abs(r.family.g - g) > 1e-12 * g {
                Issue.record("\(list.name) \(r.family.name): g \(r.family.g), pinned \(g)")
            }
        }
    }

    /// Compares one decision with its pin; returns a message when they differ.
    private func mismatch(_ plant: PinPlant, _ list: PinList, _ rings: [KinematicFamily]) -> String? {
        let k = plant.k ?? Index.defaultCoverage
        let result: Result<Index.Assignment, Index.Refusal>
        do throws(Index.Refusal) {
            result = .success(try Index.index(q: plant.q, sigma: plant.sg, in: rings, kMax: list.kmax, coverage: k))
        } catch { result = .failure(error) }
        func names(_ ids: [Int]?) -> [String] { (ids ?? []).map { list.names[$0] } }
        switch (plant.v, result) {
        case ("accepted", .success(let a)):
            guard let b = plant.b, let z = plant.z, let d = plant.d else { return "pin lacks the accepted fields" }
            if a.best.family.name != list.names[b] { return "best \(a.best.family.name), pinned \(list.names[b])" }
            if !within(a.best.z, z, 1e-10 * max(1, z)) { return "z \(a.best.z), pinned \(z)" }
            if !within(a.best.deviation, d, 1e-14) { return "delta \(a.best.deviation), pinned \(d)" }
            if let n = plant.n {
                guard let next = a.next, let margin = a.margin, let m = plant.m else { return "no next, pinned \(list.names[n])" }
                if next.family.name != list.names[n] { return "next \(next.family.name), pinned \(list.names[n])" }
                if !within(margin, m, 1e-9 * max(1, m)) { return "margin \(margin), pinned \(m)" }
            } else if a.next != nil || a.margin != nil {
                return "a next where the pin has none"
            }
            return nil
        case ("ambiguous", .failure(.ambiguous(let c))):
            let got = c.map(\.family.name), want = names(plant.p)
            return got == want ? nil : "payload \(got), pinned \(want)"
        case ("beyondKMax", .failure(.beyondKMax(let reach, let kMax))):
            guard let r = plant.reach else { return "pin lacks the reach" }
            if kMax != list.kmax { return "k_max \(kMax)" }
            return within(reach, r, 1e-9 * max(1, r)) ? nil : "reach \(reach), pinned \(r)"
        case ("noneWithinTolerance", .failure(.noneWithinTolerance(let nearest))):
            guard let nr = plant.nr, let nz = plant.nz, let nd = plant.nd else {
                return nearest == nil ? nil : "a nearest \(nearest!.family.name) where the pin has none"
            }
            guard let nearest else { return "no nearest, pinned \(list.names[nr])" }
            if nearest.family.name != list.names[nr] { return "nearest \(nearest.family.name), pinned \(list.names[nr])" }
            if !within(nearest.z, nz, 1e-10 * max(1, nz)) { return "nearest z \(nearest.z), pinned \(nz)" }
            return within(nearest.deviation, nd, 1e-14) ? nil : "nearest delta \(nearest.deviation), pinned \(nd)"
        case ("invalidInput", .failure(.invalidInput)):
            return nil
        default:
            switch result {
            case .success: return "accepted, pinned \(plant.v)"
            case .failure(let e): return "\(e.name), pinned \(plant.v)"
            }
        }
    }

    /// T1: every family of the list, u in {1e-4, 1e-3, 3e-3, 1e-2, 3e-2}, s in {+-2.9, +-2, +-1, +-0.5, 0}: q = g (1 + s u),
    /// sigma = u q, k = 3. Verdict, payload family names, delta, z, margin, reach equal to the NumPy decision. The closed
    /// forms P1 are checked on the Swift answers: an own-family acceptance has z = |s|/(1 + s u), delta = -s u/(1 + s u).
    /// Mutants: z <= k to z < k, the margin clause strict to non-strict, reach without the next family.
    @Test(arguments: ["Au", "Al", "Si", "Mg", "Mg-ideal", "Weak-cut", "beta''"])
    func testT1PlantsMatchNumPy(_ name: String) throws {
        let pins = try IndexPinFixtures.pins()
        let index = try #require(pins.t1Lists.firstIndex(of: name))
        let list = try #require(pins.lists.first { $0.name == name })
        let rings = try IndexPinFixtures.rings(list)
        var checked = 0, bad: [String] = [], p1 = 0
        for plant in pins.t1 where plant.c == index {
            checked += 1
            if let msg = mismatch(plant, list, rings) {
                bad.append("\(name) t \(plant.t ?? -1) u \(plant.u ?? 0) s \(plant.s ?? 0): \(msg)")
                continue
            }
            if plant.v == "accepted", plant.b == plant.t, let u = plant.u, let s = plant.s, let z = plant.z, let d = plant.d {
                p1 += 1
                if !within(z, abs(s) / (1 + s * u), 1e-9 * max(1, z)) || !within(d, -s * u / (1 + s * u), 1e-14) {
                    bad.append("P1 \(name) t \(plant.t ?? -1) u \(u) s \(s)")
                }
            }
        }
        #expect(checked == list.g.count * 45, "\(name): \(checked) plants")
        #expect(p1 > 0)
        #expect(bad.isEmpty, "\(bad.count) of \(checked) differ, first: \(bad.prefix(5))")
    }

    /// the reference implementation's predicted T1 counts (accepted as own, accepted as another, ambiguous, beyondKMax, none) per u over its six
    /// lists, and beta'' accepted per u (2584, 508, 184, 58, 25), reproduced from the Swift decision on the kit's rings.
    @Test func testT1CountsReproduceReferencePrediction() throws {
        let pins = try IndexPinFixtures.pins()
        let predicted: [Double: [Int]] = [1e-4: [711, 0, 0, 36, 0], 1e-3: [691, 0, 20, 36, 0], 3e-3: [647, 0, 64, 36, 0],
                                           1e-2: [408, 0, 309, 30, 0], 3e-2: [131, 19, 572, 11, 14]]
        let beta: [Double: Int] = [1e-4: 2584, 1e-3: 508, 3e-3: 184, 1e-2: 58, 3e-2: 25]
        var six: [Double: [Int]] = [:], betaAccepted: [Double: Int] = [:]
        for (index, name) in pins.t1Lists.enumerated() {
            let list = try #require(pins.lists.first { $0.name == name })
            let rings = try IndexPinFixtures.rings(list)
            for plant in pins.t1 where plant.c == index {
                let u = try #require(plant.u)
                var verdictIndex: Int
                do throws(Index.Refusal) {
                    let a = try Index.index(q: plant.q, sigma: plant.sg, in: rings, kMax: list.kmax)
                    verdictIndex = a.best.family.name == list.names[plant.t ?? -1] ? 0 : 1
                    if name == "beta''" { betaAccepted[u, default: 0] += 1 }
                } catch {
                    switch error {
                    case .ambiguous: verdictIndex = 2
                    case .beyondKMax: verdictIndex = 3
                    case .noneWithinTolerance: verdictIndex = 4
                    case .invalidInput: verdictIndex = 5
                    }
                }
                if name != "beta''", verdictIndex < 5 { six[u, default: [0, 0, 0, 0, 0]][verdictIndex] += 1 }
            }
        }
        for (u, want) in predicted { #expect(six[u] == want, "u \(u): \(String(describing: six[u])), reference predicted \(want)") }
        for (u, want) in beta { #expect(betaAccepted[u] == want, "beta'' u \(u): \(String(describing: betaAccepted[u])), predicted \(want)") }
    }

    /// T4: the refusal plants on the kit's rings: invalid numbers (q, sigma, k each 0, -1, NaN, +-inf), beyondKMax (Mg
    /// {20-25}, Si {620}, Au at 1.25), the one-family and the empty list, the forbidden positions, the ambiguous pairs (Mg,
    /// Au {331}, the margin plant, beta'' {3 1 5}/{7 3 0}) and the accepted beta'' {2 0 0}: each carries the verdict it was
    /// planted for (`expect`) and the Swift decision equals the NumPy one. Mutants: any two checks swapped.
    @Test func testT4RefusalPlantsMatchNumPy() throws {
        let pins = try IndexPinFixtures.pins()
        var bad: [String] = []
        var verdicts: [String: Int] = [:]
        for plant in pins.t4 {
            let list = try #require(pins.lists.first { $0.name == plant.l })
            let rings = try IndexPinFixtures.rings(list)
            verdicts[plant.v, default: 0] += 1
            if plant.expect != plant.v { bad.append("\(plant.label ?? "?"): pinned \(plant.v), planted for \(plant.expect ?? "?")") }
            if let msg = mismatch(plant, list, rings) { bad.append("\(plant.label ?? "?"): \(msg)") }
        }
        #expect(bad.isEmpty, "\(bad.count) differ, first: \(bad.prefix(5))")
        // every verdict is planted at least once
        for v in ["accepted", "ambiguous", "beyondKMax", "noneWithinTolerance", "invalidInput"] {
            #expect((verdicts[v] ?? 0) > 0, "no T4 plant for \(v)")
        }
        #expect(verdicts["invalidInput"] == 15)
    }
}
