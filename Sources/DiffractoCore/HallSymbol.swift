//
//  HallSymbol.swift
//  Role: the operations of a space group from its Hall symbol (Hall 1981, Acta Cryst. A37, 517; the explicit-origin
//        notation `_space_group_name_Hall` carries in a CIF). The symbol states the lattice centring, the generators
//        with their axes and translations, and an origin shift, so it fixes the group in the file's own setting with
//        no table: the 530 settings spglib 2.8.0 lists are each read to exactly its operation set
//        (`spglib.get_symmetry_from_database`, pinned by `tools/pins/cif-symmetry/hall-pins.py`).
//        Arithmetic is exact: translations are kept in twelfths.
//

import Foundation

public enum HallSymbolError: Error, Equatable, Sendable, CustomStringConvertible {
    /// The text is not a Hall symbol this reader understands; `reason` names the part that failed.
    case unreadable(symbol: String, reason: String)

    public var description: String {
        switch self {
        case .unreadable(let symbol, let reason):
            return "Could not read the Hall symbol \"\(symbol)\": \(reason)."
        }
    }
}

/// The space group a Hall symbol states, in the setting it states.
public enum HallSymbol {

    /// Every operation of the group `symbol` states, centring translations included, each translation in [0, 1). The
    /// first operation is the identity; the order is the generation order (deterministic, not sorted).
    ///
    /// The grammar is Hall's: an optional `-` (a centre of symmetry at the origin) and the lattice letter (P A B C I R S T
    /// F), then one to four matrix symbols `[-]N[axis][translations]` (N in 1 2 3 4 6; axis x y z, `'` or `"` for a
    /// two-fold on a face diagonal, `*` for the body diagonal; translations a b c n u v w d and screw digits), then an
    /// optional origin shift `(i j k)` in twelfths. An axis left out takes Hall's default: z first, x for a two-fold
    /// after a 2 or 4, `'` for a two-fold after a 3 or 6, `*` for a three-fold in third place. Lattice letters are read
    /// in either case; everything else is case-sensitive, as Hall wrote it.
    public static func operations(_ symbol: String) throws(HallSymbolError) -> [SymmetryOperation] {
        try group(symbol).map { op in
            SymmetryOperation(rotation: op.rotation,
                              translation: SIMD3(Double(op.twelfths.x), Double(op.twelfths.y), Double(op.twelfths.z)) / 12)
        }
    }

    // MARK: - Exact operations

    /// R p + t with t in twelfths, reduced mod 12.
    struct ExactOperation: Hashable {
        let rotation: IntegerMatrix
        let twelfths: SIMD3<Int>

        init(rotation: IntegerMatrix, twelfths: SIMD3<Int>) {
            self.rotation = rotation
            self.twelfths = SIMD3(((twelfths.x % 12) + 12) % 12, ((twelfths.y % 12) + 12) % 12,
                                  ((twelfths.z % 12) + 12) % 12)
        }

        static let identity = ExactOperation(rotation: .identity, twelfths: .zero)

        /// self after other: R₁R₂, R₁t₂ + t₁.
        func composed(with other: ExactOperation) -> ExactOperation {
            ExactOperation(rotation: rotation.multiplied(by: other.rotation),
                           twelfths: HallSymbol.apply(rotation, other.twelfths) &+ twelfths)
        }
    }

    static func apply(_ m: IntegerMatrix, _ v: SIMD3<Int>) -> SIMD3<Int> {
        SIMD3((m.rows[0] &* v).wrappedSum(), (m.rows[1] &* v).wrappedSum(), (m.rows[2] &* v).wrappedSum())
    }

    // MARK: - Tables (Hall 1981, tables 2 to 4)

    /// Lattice translations per letter, in twelfths, beyond the origin.
    private static let centring: [Character: [SIMD3<Int>]] = [
        "P": [], "A": [SIMD3(0, 6, 6)], "B": [SIMD3(6, 0, 6)], "C": [SIMD3(6, 6, 0)], "I": [SIMD3(6, 6, 6)],
        "R": [SIMD3(8, 4, 4), SIMD3(4, 8, 8)], "S": [SIMD3(4, 4, 8), SIMD3(8, 8, 4)],
        "T": [SIMD3(4, 8, 4), SIMD3(8, 4, 8)],
        "F": [SIMD3(0, 6, 6), SIMD3(6, 0, 6), SIMD3(6, 6, 0)],
    ]

    /// Translation symbols, in twelfths.
    private static let translationSymbols: [Character: SIMD3<Int>] = [
        "a": SIMD3(6, 0, 0), "b": SIMD3(0, 6, 0), "c": SIMD3(0, 0, 6), "n": SIMD3(6, 6, 6),
        "u": SIMD3(3, 0, 0), "v": SIMD3(0, 3, 0), "w": SIMD3(0, 0, 3), "d": SIMD3(3, 3, 3),
    ]

    /// The proper rotation of order `n` about a principal axis.
    private static func principal(_ n: Int, axis: Character) -> IntegerMatrix? {
        let z: IntegerMatrix
        switch n {
        case 1: return .identity
        case 2: z = IntegerMatrix([-1, 0, 0], [0, -1, 0], [0, 0, 1])
        case 3: z = IntegerMatrix([0, -1, 0], [1, -1, 0], [0, 0, 1])
        case 4: z = IntegerMatrix([0, -1, 0], [1, 0, 0], [0, 0, 1])
        case 6: z = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        default: return nil
        }
        // x and y by cyclic permutation of the axes: (x, y, z) -> (y, z, x) takes the z-axis matrix to the x-axis one.
        func permuted(_ m: IntegerMatrix, _ p: [Int]) -> IntegerMatrix {
            IntegerMatrix(rows: (0..<3).map { r in SIMD3((0..<3).map { c in m[p[r], p[c]] }) })
        }
        switch axis {
        case "z": return z
        case "x": return permuted(z, [2, 0, 1])
        case "y": return permuted(z, [1, 2, 0])
        default: return nil
        }
    }

    /// The two-folds on face diagonals (`'` along b - a style, `"` along a + b style), relative to the preceding axis.
    private static func diagonalTwoFold(_ kind: Character, preceding: Character) -> IntegerMatrix? {
        switch (kind, preceding) {
        case ("'", "z"): return IntegerMatrix([0, -1, 0], [-1, 0, 0], [0, 0, -1])
        case ("\"", "z"): return IntegerMatrix([0, 1, 0], [1, 0, 0], [0, 0, -1])
        case ("'", "x"): return IntegerMatrix([-1, 0, 0], [0, 0, -1], [0, -1, 0])
        case ("\"", "x"): return IntegerMatrix([-1, 0, 0], [0, 0, 1], [0, 1, 0])
        case ("'", "y"): return IntegerMatrix([0, 0, -1], [0, -1, 0], [-1, 0, 0])
        case ("\"", "y"): return IntegerMatrix([0, 0, 1], [0, -1, 0], [1, 0, 0])
        default: return nil
        }
    }

    /// The three-fold along the body diagonal.
    private static let bodyDiagonalThreeFold = IntegerMatrix([0, 0, 1], [1, 0, 0], [0, 1, 0])

    // MARK: - Reading

    /// The group as exact operations, the identity first.
    static func group(_ symbol: String) throws(HallSymbolError) -> [ExactOperation] {
        func fail(_ reason: String) -> HallSymbolError { .unreadable(symbol: symbol, reason: reason) }

        var text = symbol.trimmingCharacters(in: .whitespaces)
        var shift = SIMD3<Int>.zero
        if let open = text.firstIndex(of: "(") {
            guard text.hasSuffix(")") else { throw fail("the origin shift is not closed by \")\"") }
            let inner = text[text.index(after: open)..<text.index(before: text.endIndex)]
            let parts = inner.split(whereSeparator: { $0.isWhitespace })
            guard parts.count == 3, let i = Int(parts[0]), let j = Int(parts[1]), let k = Int(parts[2]) else {
                throw fail("the origin shift \"(\(inner))\" is not three whole numbers of twelfths")
            }
            shift = SIMD3(i % 12, j % 12, k % 12)
            text = String(text[text.startIndex..<open]).trimmingCharacters(in: .whitespaces)
        }

        var tokens = text.split(whereSeparator: { $0.isWhitespace }).map(String.init)
        guard !tokens.isEmpty else { throw fail("it is empty") }
        // The lattice symbol, written alone or run into the first matrix ("-P2ybc").
        var latticeToken = tokens.removeFirst()
        let centric = latticeToken.hasPrefix("-")
        if centric { latticeToken.removeFirst() }
        guard let first = latticeToken.first, first.isASCII, let letter = first.uppercased().first,
              let centringVectors = centring[letter]
        else { throw fail("it does not start with a lattice letter (P A B C I R S T F)") }
        let attached = String(latticeToken.dropFirst())
        if !attached.isEmpty { tokens.insert(attached, at: 0) }
        guard !tokens.isEmpty else { throw fail("it has no matrix symbol") }
        guard tokens.count <= 4 else { throw fail("it has more than four matrix symbols") }

        var generators: [ExactOperation] = centringVectors.map { ExactOperation(rotation: .identity, twelfths: $0) }
        if centric { generators.append(ExactOperation(rotation: .inversion, twelfths: .zero)) }

        var previousOrder = 0
        var previousAxis: Character = "z"
        for (index, token) in tokens.enumerated() {
            var chars = Array(token)
            var improper = false
            if chars.first == "-" { improper = true; chars.removeFirst() }
            guard let first = chars.first, let n = first.wholeNumberValue, [1, 2, 3, 4, 6].contains(n) else {
                throw fail("\"\(token)\" does not start with a rotation order (1 2 3 4 6)")
            }
            chars.removeFirst()
            // Hall writes a 1 after the first place only as -1 or with a translation; a bare "1" there is a
            // Hermann–Mauguin name typed into the Hall tag ("P 3 2 1" is P321, which Hall writes "P 3 2\"").
            if n == 1, index > 0, !improper, chars.isEmpty {
                throw fail("\"\(token)\" is a bare 1 after the first place (a Hermann–Mauguin symbol, not a Hall symbol?)")
            }
            var axis: Character?
            if let c = chars.first, "xyz'\"*".contains(c) { axis = c; chars.removeFirst() }
            if axis == nil {
                switch (index, n) {
                case (0, _): axis = "z"
                case (1, 2) where previousOrder == 2 || previousOrder == 4: axis = "x"
                case (1, 2) where previousOrder == 3 || previousOrder == 6: axis = "'"
                case (2, 3): axis = "*"
                case (_, 1): axis = "z"
                default: throw fail("\"\(token)\" leaves out an axis Hall gives no default for")
                }
            }
            guard let axis else { throw fail("\"\(token)\" has no axis") }

            let proper: IntegerMatrix?
            switch axis {
            case "x", "y", "z": proper = principal(n, axis: axis)
            case "'", "\"": proper = n == 2 ? diagonalTwoFold(axis, preceding: previousAxis) : nil
            default: proper = n == 3 ? bodyDiagonalThreeFold : nil
            }
            guard let proper else { throw fail("\"\(token)\" names an axis its rotation order cannot have") }

            var translation = SIMD3<Int>.zero
            for c in chars {
                if let t = translationSymbols[c] {
                    translation &+= t
                } else if let k = c.wholeNumberValue, k >= 1, k < n, "xyz".contains(axis) {
                    let step = 12 * k / n
                    switch axis {
                    case "x": translation.x += step
                    case "y": translation.y += step
                    default: translation.z += step
                    }
                } else {
                    throw fail("\"\(token)\" has a translation symbol \"\(c)\" that is not a b c n u v w d or a screw digit below \(n)")
                }
            }
            generators.append(ExactOperation(rotation: improper ? proper.negated : proper, twelfths: translation))
            if "xyz".contains(axis) { previousAxis = axis }
            previousOrder = n
        }

        // The closure; no space group has more than 192 operations in any Hall setting.
        var elements: [ExactOperation] = [.identity]
        var seen: Set<ExactOperation> = [.identity]
        var frontier = 0
        while frontier < elements.count {
            let element = elements[frontier]
            frontier += 1
            for generator in generators {
                let product = generator.composed(with: element)
                if seen.insert(product).inserted {
                    elements.append(product)
                    guard elements.count <= 192 else { throw fail("its generators do not close into a space group") }
                }
            }
        }

        // The origin shift: the operation x -> R x + t seen from an origin moved by -v is x -> R x + t + (I - R) v.
        guard shift != .zero else { return elements }
        return elements.map { op in
            ExactOperation(rotation: op.rotation, twelfths: op.twelfths &+ shift &- apply(op.rotation, shift))
        }
    }
}
