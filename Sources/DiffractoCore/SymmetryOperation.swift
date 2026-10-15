//
//  SymmetryOperation.swift
//  Role: one affine operation on fractional coordinates, p' = R p + t, with R an integer matrix and t a translation;
//        parsed from and printed as the CIF `x,y,z` grammar. The parser is reference commit 7448401b's
//        `CIFImport.parseSymmetryOperation` (registration section 3.3), with one addition: a rotation coefficient that
//        is not an integer is refused, since no space-group operation has one in any setting.
//

import Foundation

public enum SymmetryOperationError: Error, Equatable, Sendable, CustomStringConvertible {
    /// The text is not three comma-separated affine expressions in x, y, z.
    case unparsable(String)
    /// A coefficient of x, y or z is not an integer.
    case nonIntegerRotation(String)
    /// The rotation part's determinant is not ±1 (`x,x,z`, `2x,y,z`): no symmetry operation is singular or a dilation.
    case notInvertible(String)

    public var description: String {
        switch self {
        case .unparsable(let raw): return "Could not parse symmetry operation \"\(raw)\" as an x,y,z expression."
        case .nonIntegerRotation(let raw):
            return "Symmetry operation \"\(raw)\" has a non-integer coefficient; the rotation part of a space-group "
                + "operation is an integer matrix in every setting."
        case .notInvertible(let raw):
            return "Symmetry operation \"\(raw)\" is not a symmetry: its rotation part has a determinant other than ±1."
        }
    }
}

/// p' = R p + t on fractional coordinates.
public struct SymmetryOperation: Sendable, Hashable, CustomStringConvertible {
    public let rotation: IntegerMatrix
    public let translation: SIMD3<Double>

    public init(rotation: IntegerMatrix, translation: SIMD3<Double> = SIMD3(0, 0, 0)) {
        self.rotation = rotation
        self.translation = translation
    }

    public static let identity = SymmetryOperation(rotation: .identity)

    /// Parses `"x,y,z"`, `"-x,-y,-z"`, `"1/2+x,1/2-y,z"`, `"x-y,x,-z+1/3"`, `"-X, Y+1/2, -Z"`, …
    public init(xyz raw: String) throws(SymmetryOperationError) {
        let components = raw.split(separator: ",", omittingEmptySubsequences: false)
            .map { $0.trimmingCharacters(in: .whitespaces) }
        guard components.count == 3, components.allSatisfy({ !$0.isEmpty }) else {
            throw .unparsable(raw)
        }
        var rows: [SIMD3<Int>] = []
        var translation = SIMD3<Double>(0, 0, 0)
        for (axis, component) in components.enumerated() {
            let (coefficients, constant) = try SymmetryOperation.parseComponent(component, raw: raw)
            var row = SIMD3<Int>(0, 0, 0)
            for i in 0..<3 {
                let rounded = coefficients[i].rounded()
                guard abs(coefficients[i] - rounded) <= 1e-9 else { throw .nonIntegerRotation(raw) }
                row[i] = Int(rounded)
            }
            rows.append(row)
            translation[axis] = constant
        }
        let rotation = IntegerMatrix(rows: rows)
        guard abs(rotation.determinant) == 1 else { throw .notInvertible(raw) }
        self.rotation = rotation
        self.translation = translation
    }

    /// Parses one component (e.g. `"1/2+x"`) into the (x, y, z) coefficients and the constant.
    private static func parseComponent(
        _ component: String, raw: String
    ) throws(SymmetryOperationError) -> (SIMD3<Double>, Double) {
        let chars = Array(component.lowercased().replacingOccurrences(of: " ", with: ""))
        guard !chars.isEmpty else { throw .unparsable(raw) }

        var coefficients = SIMD3<Double>(0, 0, 0)
        var constant = 0.0
        var i = 0
        while i < chars.count {
            var sign = 1.0
            if chars[i] == "+" { i += 1 } else if chars[i] == "-" { sign = -1.0; i += 1 }
            guard i < chars.count else { throw .unparsable(raw) }

            let numStart = i
            while i < chars.count, chars[i].isNumber || chars[i] == "/" || chars[i] == "." {
                i += 1
            }
            let numericString = String(chars[numStart..<i])
            let hasNumeric = !numericString.isEmpty
            let magnitude = hasNumeric ? try parseFraction(numericString, raw: raw) : 1.0

            if i < chars.count, "xyz".contains(chars[i]) {
                switch chars[i] {
                case "x": coefficients.x += sign * magnitude
                case "y": coefficients.y += sign * magnitude
                default: coefficients.z += sign * magnitude
                }
                i += 1
                // A term ends at a sign or at the end: `xy`, `x2` and `x.5` are refused, not read as x + y, x + 2,
                // x + 0.5 (reference commit 7d6fdbd1's rule, its refuter RF5, 2026-10-15).
                if i < chars.count, chars[i] != "+", chars[i] != "-" { throw .unparsable(raw) }
            } else if hasNumeric {
                constant += sign * magnitude
            } else {
                throw .unparsable(raw)
            }
        }
        return (coefficients, constant)
    }

    /// `"1"`, `"0.5"`, `"1/2"`; an empty part (`"/2"`), a zero denominator, a non-finite value or a magnitude beyond
    /// 1000 (a garbled `1000000000000000000000x`) is unparsable, never a trap.
    private static func parseFraction(_ text: String, raw: String) throws(SymmetryOperationError) -> Double {
        let parts = text.split(separator: "/", omittingEmptySubsequences: false)
        var value: Double?
        if parts.count == 1 { value = Double(parts[0]) }
        if parts.count == 2, let numerator = Double(parts[0]), let denominator = Double(parts[1]), denominator != 0 {
            value = numerator / denominator
        }
        guard let value, value.isFinite, abs(value) <= 1000 else { throw .unparsable(raw) }
        return value
    }

    // MARK: Behaviour

    /// R p + t.
    public func apply(_ p: SIMD3<Double>) -> SIMD3<Double> {
        rotation.apply(p) + translation
    }

    /// The operation p ↦ self(other(p)): rotation R₁R₂, translation R₁t₂ + t₁.
    public func composed(with other: SymmetryOperation) -> SymmetryOperation {
        SymmetryOperation(rotation: rotation.multiplied(by: other.rotation),
                          translation: rotation.apply(other.translation) + translation)
    }

    /// The translation with each component in [0, 1).
    public var translationMod1: SIMD3<Double> {
        wrapFractional(translation)
    }

    /// True for `x,y,z` and its lattice-translated spellings (`x+1,y,z`): the operation expands nothing.
    public var isIdentity: Bool {
        rotation.isIdentity && SymmetryOperation.isWholeNumber(translation)
    }

    /// Same rotation, translations equal mod 1 within `translationTolerance`.
    public func isEquivalent(to other: SymmetryOperation, translationTolerance: Double = 1e-9) -> Bool {
        guard rotation == other.rotation else { return false }
        return fractionalAxisDelta(translation.x, other.translation.x) <= translationTolerance
            && fractionalAxisDelta(translation.y, other.translation.y) <= translationTolerance
            && fractionalAxisDelta(translation.z, other.translation.z) <= translationTolerance
    }

    /// The number of distinct operations: rotation exact, translation mod 1 to 1e-3 (the finest genuine translation
    /// is 1/12 apart), so duplicated rows cannot pad an incomplete list. A garbled coefficient (1e21) is clamped,
    /// not trapped, and counts as its own operation (the reference implementations' rule).
    public static func distinctCount(_ ops: [SymmetryOperation]) -> Int {
        func quantised(_ value: Double, scale: Double) -> Int {
            guard value.isFinite else { return Int.max }
            return Int((max(-1e9, min(1e9, value)) * scale).rounded())
        }
        var seen = Set<[Int]>()
        for op in ops {
            var key: [Int] = []
            for row in op.rotation.rows { key.append(contentsOf: [row.x, row.y, row.z]) }
            for t in [op.translation.x, op.translation.y, op.translation.z] {
                key.append(quantised(t - t.rounded(.down), scale: 1e3) % 1000)
            }
            seen.insert(key)
        }
        return seen.count
    }

    static func isWholeNumber(_ v: SIMD3<Double>) -> Bool {
        func whole(_ x: Double) -> Bool { abs(x - x.rounded()) <= 1e-9 }
        return whole(v.x) && whole(v.y) && whole(v.z)
    }

    // MARK: Printing

    /// The CIF form, `x,y,z`, `-x,y+1/2,-z`: variables first, then the translation as a fraction of 2, 3, 4, 6 or 12
    /// when it is one within 1e-9, else as a decimal.
    public var xyzString: String {
        let names = ["x", "y", "z"]
        var parts: [String] = []
        for axis in 0..<3 {
            var s = ""
            let row = rotation.rows[axis]
            for (i, coefficient) in [row.x, row.y, row.z].enumerated() where coefficient != 0 {
                let magnitude = abs(coefficient) == 1 ? "" : "\(abs(coefficient))"
                let sign = coefficient < 0 ? "-" : (s.isEmpty ? "" : "+")
                s += sign + magnitude + names[i]
            }
            let t = translation[axis]
            if abs(t) > 1e-9 {
                let fraction = SymmetryOperation.fractionText(t)
                s += (t < 0 ? "-" : (s.isEmpty ? "" : "+")) + fraction
            }
            parts.append(s.isEmpty ? "0" : s)
        }
        return parts.joined(separator: ",")
    }

    private static func fractionText(_ t: Double) -> String {
        let magnitude = abs(t)
        for denominator in [1, 2, 3, 4, 6, 12] {
            let scaled = magnitude * Double(denominator)
            if abs(scaled - scaled.rounded()) <= 1e-9 {
                let numerator = Int(scaled.rounded())
                return denominator == 1 ? "\(numerator)" : "\(numerator)/\(denominator)"
            }
        }
        return String(format: "%.6g", magnitude)
    }

    public var description: String { xyzString }
}
