//
//  CellStandardisation.swift
//  Role: the primitive-to-conventional cell conversion of a centred lattice by the fixed ITA change of basis for its
//        centring letter (F, I, C, A, B, R obverse), with the guards that refuse a conversion that does not give a
//        conventional cell of the stated crystal system. Ported from mac4DSTEM `MaterialsProjectImport.standardise`
//        (the numbers only: the JSON types and the fetch stay in the app). This is the app's route, pinned to the app's
//        output; the reference route (pymatgen SpacegroupAnalyzer) is `CrystalStructure.detectingSymmetry(...)
//        .conventionalSetting`, and both are tested against pymatgen's standard cells (docs/parity.md).
//

import Foundation
import simd

/// Conversion of a served primitive cell of a centred lattice to its conventional cell.
public enum CellStandardisation {
    /// Whether the cell had to be converted.
    public enum Outcome: Sendable, Equatable {
        /// The cell already was conventional (a `P` centring, or a centred cell that reads conventional), or no
        /// symmetry was stated.
        case asServed
        /// The cell was a primitive cell for this centring letter and was converted.
        case conventional(centring: Character)
    }

    /// A standardised cell: lattice rows (Å), sites and what was done.
    public struct Result: Sendable, Equatable {
        /// The lattice rows (Å).
        public let latticeRows: [SIMD3<Double>]
        /// The sites in the conventional cell (fractional, in [0, 1)).
        public let sites: [AtomSite]
        /// Whether a conversion happened.
        public let outcome: Outcome
    }

    /// The conversion was refused; the message says which guard failed.
    public struct Failure: Error, Sendable, Equatable, CustomStringConvertible {
        /// Which guard failed, in words.
        public let message: String
        public var description: String { message }
    }

    /// Fractional-distance floor between two sites of one element after the conversion.
    public static let minimumSiteSeparation = 1e-3

    /// Converts `latticeRows` and `sites` to the conventional cell of the centring that starts `symbol` (for example
    /// `"Fm-3m"`), checking the result against `crystalSystem` (lower-case: "cubic", "tetragonal", "orthorhombic",
    /// "hexagonal", "trigonal", "monoclinic", "triclinic").
    ///
    /// A missing symbol or a `P` centring returns the input. Otherwise a cell that already reads conventional for the
    /// system (lengths to 1e-3 relative, angles to 1e-2°) is returned as served; else the fixed matrix is applied and
    /// the result must pass: the cell reads conventional for the system, the volume is the centring's multiplicity times
    /// the served volume (1e-6 relative), the site count and each element's count scale by the multiplicity, and no two
    /// same-element sites are closer than ``minimumSiteSeparation`` in fractional distance.
    public static func standardise(
        latticeRows: [SIMD3<Double>], sites: [AtomSite], symbol: String?, crystalSystem: String?
    ) throws -> Result {
        guard latticeRows.count == 3 else { throw Failure(message: "the lattice needs three rows.") }
        guard !sites.isEmpty else { throw Failure(message: "there are no sites to standardise.") }
        guard let centring = symbol?.first else {
            return Result(latticeRows: latticeRows, sites: sites, outcome: .asServed)
        }
        if centring == "P" { return Result(latticeRows: latticeRows, sites: sites, outcome: .asServed) }
        guard let system = crystalSystem?.lowercased() else {
            throw Failure(message: "no crystal system was given to verify a \(centring)-centred cell against.")
        }
        let served = parameters(latticeRows)
        if reads(conventional: system, served) {
            return Result(latticeRows: latticeRows, sites: sites, outcome: .asServed)
        }
        guard let transform = transform(for: centring) else {
            throw Failure(message: "unrecognised centring letter '\(centring)'.")
        }
        guard let inverse = invert(transform.matrix) else {
            throw Failure(message: "the \(centring)-centring matrix is not invertible.")
        }
        let rows = transformRows(latticeRows, by: transform.matrix)
        let new = parameters(rows)
        guard reads(conventional: system, new) else {
            throw Failure(message: "converting as \(centring)-centred did not give a conventional \(system) cell "
                + "(a=\(new.a), b=\(new.b), c=\(new.c), α=\(new.alpha)°, β=\(new.beta)°, γ=\(new.gamma)°).")
        }
        let servedVolume = abs(determinant(latticeRows))
        guard servedVolume > 0, servedVolume.isFinite else {
            throw Failure(message: "the served cell has no finite positive volume.")
        }
        let ratio = abs(determinant(rows)) / servedVolume
        let multiplicity = Double(transform.translations.count)
        guard abs(ratio - multiplicity) <= 1e-6 * multiplicity else {
            throw Failure(message: "the converted volume is \(ratio)x the served one, not \(multiplicity)x.")
        }
        var out: [AtomSite] = []
        out.reserveCapacity(sites.count * transform.translations.count)
        for site in sites {
            let base = rowTimes(site.fractional, inverse)
            for t in transform.translations {
                var s = site
                s.fractional = wrap(base + t)
                out.append(s)
            }
        }
        var servedCounts: [Int: Int] = [:], newCounts: [Int: Int] = [:]
        for s in sites { servedCounts[s.z, default: 0] += 1 }
        for s in out { newCounts[s.z, default: 0] += 1 }
        for (z, n) in servedCounts where newCounts[z] != n * transform.translations.count {
            throw Failure(message: "the site count of Z=\(z) did not scale by \(transform.translations.count).")
        }
        for i in 0..<out.count {
            for j in (i + 1)..<out.count where out[i].z == out[j].z {
                let d = out[i].fractional - out[j].fractional
                let wrapped = SIMD3(axisDelta(d.x), axisDelta(d.y), axisDelta(d.z))
                if (wrapped * wrapped).sum().squareRoot() < minimumSiteSeparation {
                    throw Failure(message: "two sites of Z=\(out[i].z) landed on each other after the conversion; "
                        + "the served sites are not this centring's asymmetric unit.")
                }
            }
        }
        return Result(latticeRows: rows, sites: out, outcome: .conventional(centring: centring))
    }

    // MARK: - Private

    private struct Transform { let matrix: [[Double]]; let translations: [SIMD3<Double>] }

    private static func transform(for letter: Character) -> Transform? {
        switch letter {
        case "F": return Transform(matrix: [[-1, 1, 1], [1, -1, 1], [1, 1, -1]],
                                   translations: [SIMD3(0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)])
        case "I": return Transform(matrix: [[0, 1, 1], [1, 0, 1], [1, 1, 0]], translations: [SIMD3(0, 0, 0), SIMD3(0.5, 0.5, 0.5)])
        case "C": return Transform(matrix: [[1, -1, 0], [1, 1, 0], [0, 0, 1]], translations: [SIMD3(0, 0, 0), SIMD3(0.5, 0.5, 0)])
        case "A": return Transform(matrix: [[1, 0, 0], [0, 1, -1], [0, 1, 1]], translations: [SIMD3(0, 0, 0), SIMD3(0, 0.5, 0.5)])
        case "B": return Transform(matrix: [[1, 0, -1], [0, 1, 0], [1, 0, 1]], translations: [SIMD3(0, 0, 0), SIMD3(0.5, 0, 0.5)])
        case "R": return Transform(matrix: [[1, -1, 0], [0, 1, -1], [1, 1, 1]],
                                   translations: [SIMD3(0, 0, 0), SIMD3(2.0 / 3, 1.0 / 3, 1.0 / 3), SIMD3(1.0 / 3, 2.0 / 3, 2.0 / 3)])
        default: return nil
        }
    }

    private struct Params { let a, b, c, alpha, beta, gamma: Double }

    private static func angle(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double {
        acos(min(1, max(-1, simd_dot(u, v) / (simd_length(u) * simd_length(v))))) * 180 / .pi
    }

    private static func parameters(_ r: [SIMD3<Double>]) -> Params {
        Params(a: simd_length(r[0]), b: simd_length(r[1]), c: simd_length(r[2]),
               alpha: angle(r[1], r[2]), beta: angle(r[0], r[2]), gamma: angle(r[0], r[1]))
    }

    private static func reads(conventional system: String, _ p: Params) -> Bool {
        func len(_ x: Double, _ y: Double) -> Bool { abs(x - y) <= 1e-3 * max(1, max(abs(x), abs(y))) }
        func ang(_ x: Double, _ t: Double) -> Bool { abs(x - t) <= 1e-2 }
        switch system {
        case "cubic": return len(p.a, p.b) && len(p.a, p.c) && ang(p.alpha, 90) && ang(p.beta, 90) && ang(p.gamma, 90)
        case "tetragonal": return len(p.a, p.b) && ang(p.alpha, 90) && ang(p.beta, 90) && ang(p.gamma, 90)
        case "orthorhombic": return ang(p.alpha, 90) && ang(p.beta, 90) && ang(p.gamma, 90)
        case "hexagonal", "trigonal": return len(p.a, p.b) && ang(p.alpha, 90) && ang(p.beta, 90) && ang(p.gamma, 120)
        case "monoclinic": return ang(p.alpha, 90) && ang(p.gamma, 90)
        case "triclinic": return true
        default: return false
        }
    }

    private static func determinant(_ m: [[Double]]) -> Double {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    private static func determinant(_ r: [SIMD3<Double>]) -> Double {
        determinant(r.map { [$0.x, $0.y, $0.z] })
    }

    private static func invert(_ m: [[Double]]) -> [[Double]]? {
        let det = determinant(m)
        guard abs(det) > 1e-12 else { return nil }
        let cof: [[Double]] = [
            [m[1][1] * m[2][2] - m[1][2] * m[2][1], -(m[1][0] * m[2][2] - m[1][2] * m[2][0]), m[1][0] * m[2][1] - m[1][1] * m[2][0]],
            [-(m[0][1] * m[2][2] - m[0][2] * m[2][1]), m[0][0] * m[2][2] - m[0][2] * m[2][0], -(m[0][0] * m[2][1] - m[0][1] * m[2][0])],
            [m[0][1] * m[1][2] - m[0][2] * m[1][1], -(m[0][0] * m[1][2] - m[0][2] * m[1][0]), m[0][0] * m[1][1] - m[0][1] * m[1][0]],
        ]
        var inv = [[Double]](repeating: [0, 0, 0], count: 3)
        for i in 0..<3 { for j in 0..<3 { inv[i][j] = cof[j][i] / det } }
        return inv
    }

    /// New row i = sum over j of m[i][j] times old row j, summed left to right.
    private static func transformRows(_ old: [SIMD3<Double>], by m: [[Double]]) -> [SIMD3<Double>] {
        (0..<3).map { i in
            let t0: SIMD3<Double> = m[i][0] * old[0]
            let t1: SIMD3<Double> = m[i][1] * old[1]
            let t2: SIMD3<Double> = m[i][2] * old[2]
            let partial: SIMD3<Double> = t0 + t1
            return partial + t2
        }
    }

    private static func rowTimes(_ v: SIMD3<Double>, _ m: [[Double]]) -> SIMD3<Double> {
        SIMD3(v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
              v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
              v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2])
    }

    private static func wrap(_ v: SIMD3<Double>) -> SIMD3<Double> {
        func one(_ x: Double) -> Double {
            var m = x.truncatingRemainder(dividingBy: 1)
            if m < 0 { m += 1 }
            if m >= 1 - 1e-9 || abs(m) < 1e-9 { m = 0 }
            return m
        }
        return SIMD3(one(v.x), one(v.y), one(v.z))
    }

    private static func axisDelta(_ d: Double) -> Double {
        let m = abs(d).truncatingRemainder(dividingBy: 1)
        return min(m, 1 - m)
    }
}
