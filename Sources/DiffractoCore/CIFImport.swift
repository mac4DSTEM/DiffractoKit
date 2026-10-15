//
//  CIFImport.swift
//  Role: the subset of CIF that defines a structure (cell, the atom-site loop, the symmetry-operation loop), read into
//        a `CrystalStructure`. Ported from reference commit 7448401b `CIFImport.swift` (registration section 3.7): the
//        tokenizer, the `(3)` uncertainty suffix, the element symbol rule, the written-precision half-step, the
//        expansion that merges images within the propagated half-step AND 0.10 Å, and the refusals that protect the
//        structure. Generalised: every crystal system is admitted; a file without operations has its symmetry
//        detected from the sites (3.4b); the two app-specific refusals (point group, scattering factor) are gone.
//        Errors are user-facing: every case names the offending tag, value or symbol so a scientist can find and fix
//        the problem in their own file.
//

import Foundation
import simd

public enum CIFImportError: Error, Equatable, Sendable, CustomStringConvertible {
    /// A required cell-parameter or atom-site coordinate tag is absent.
    case missingTag(String)
    /// A numeric value could not be parsed, even after stripping a standard-uncertainty suffix like `(3)`.
    case unparsableNumber(tag: String, value: String)
    /// No `_atom_site_*` loop was found, or it produced zero rows.
    case noAtomSites
    /// An atom-site label or type symbol did not resolve to a chemical element.
    case unknownElementSymbol(String, atomSiteLabel: String)
    /// A `_symmetry_equiv_pos_as_xyz` / `_space_group_symop_operation_xyz` entry is not an affine x,y,z expression with
    /// integer coefficients.
    case unsupportedSymmetryOperation(String)
    /// The file declares a space group other than P1 but carries no symmetry operations to expand the asymmetric unit
    /// with (or only the identity). Importing the listed sites as the complete cell would silently build a wrong
    /// crystal (rock salt read this way becomes the CsCl structure).
    case missingSymmetryOperations(spaceGroup: String)
    /// A `loop_`'s value count is not a whole multiple of its tag count: the file is truncated or malformed, and
    /// continuing would silently drop the partial trailing row.
    case truncatedLoop(firstTag: String, columns: Int, values: Int)
    /// The assembled structure failed `CrystalStructure.validationIssues`.
    case invalidStructure([String])
    /// More than one `data_` block defines a structure (a cell or an atom-site loop).
    case multipleStructures(blockNames: [String])
    /// The file names its space group by IT number but lists fewer distinct symmetry operations than that group has
    /// in the conventional setting; expanding with a partial list builds a smaller cell content.
    case incompleteSymmetryOperations(spaceGroupNumber: Int, declaredName: String?, listed: Int, required: Int)
    /// A `data_` block that defines no structure carries a symmetry-operation loop while another block defines one.
    case symmetryOperationsOutsideStructure(blockNames: [String])

    public var description: String {
        switch self {
        case .missingTag(let tag):
            return "CIF is missing the required tag \(tag)."
        case .unparsableNumber(let tag, let value):
            return "Could not parse a number from \(tag) = \"\(value)\"."
        case .noAtomSites:
            return "CIF has no _atom_site_* loop (or it has no rows) — there is no atomic basis to import."
        case .unknownElementSymbol(let symbol, let label):
            return "Atom site \"\(label)\" has an unrecognized element symbol \"\(symbol)\"."
        case .unsupportedSymmetryOperation(let op):
            return "Could not parse symmetry operation \"\(op)\" as an x,y,z expression with integer coefficients."
        case .missingSymmetryOperations(let group):
            return "CIF declares space group \(group) but lists no symmetry operations beyond the identity, so its atom "
                + "sites are only the asymmetric unit. Importing them as the complete cell would build a wrong crystal, "
                + "so the file is refused. Re-export it with its symmetry-operation loop (_symmetry_equiv_pos_as_xyz or "
                + "_space_group_symop_operation_xyz — VESTA and pymatgen both write one), or supply a P1 file that lists "
                + "every atom in the cell."
        case .truncatedLoop(let firstTag, let columns, let values):
            return "The loop starting at \(firstTag) has \(columns) column(s) but \(values) value(s) — not a whole number "
                + "of rows, so the file is truncated or malformed. It is refused rather than silently dropping the "
                + "partial trailing row."
        case .invalidStructure(let issues):
            return "Imported crystal structure failed validation: \(issues.joined(separator: "; "))."
        case .multipleStructures(let names):
            return "This CIF holds \(names.count) structures (\(names.joined(separator: ", "))); import one at a time. "
                + "Split the file so each data_ block is its own CIF — reading them together would merge their cells "
                + "and atoms into one crystal that is none of them."
        case .incompleteSymmetryOperations(let number, let name, let listed, let required):
            let group = name.map { "No. \(number) (\($0))" } ?? "No. \(number)"
            return "The CIF lists \(listed) distinct symmetry operations but space group \(group) has \(required) in its "
                + "conventional setting — the list is incomplete (or the space-group number does not match the "
                + "operations), so expanding the atom sites with it would build a smaller structure than the file "
                + "describes (rock salt without its F-centring reads as CsCl). Re-export the CIF with its full "
                + "symmetry-operation loop, or supply a P1 file that lists every atom in the cell."
        case .symmetryOperationsOutsideStructure(let names):
            return "The CIF has a symmetry-operation loop in a data_ block that holds no structure "
                + "(\(names.joined(separator: ", "))). Reading the file as one crystal would apply those operations to "
                + "another block's atoms. Move the loop into the structure's own data_ block, or split the file."
        }
    }
}

/// `localizedDescription` (what an app shows) is the same text as `description` (registration 0.0.4 section 5).
extension CIFImportError: LocalizedError {
    public var errorDescription: String? { description }
}

/// What a CIF read yields: the structure and what the file said about itself.
public struct CIFDocument: Sendable, Equatable {
    public let structure: CrystalStructure
    /// The first `data_` block's name; nil when it is empty.
    public let dataBlockName: String?
    /// `CIFImport.displayName`: the block name unless it is empty or a structure editor's placeholder.
    public let displayName: String
    /// The Hermann–Mauguin symbol the file states, as written; nil when none.
    public let declaredSpaceGroupSymbol: String?
    /// The IT number the file states; nil when none.
    public let declaredSpaceGroupNumber: Int?
    /// The coarsest rounding half-step over every written fractional coordinate (0 when all are exact).
    public let coarsestHalfStep: Double
    /// Whether the file carried an expanding symmetry-operation loop (`structure.fileOperations`); the structure's
    /// operations are detected from the expanded cell either way.
    public let hadSymmetryOperations: Bool
}

public enum CIFImport {

    /// Reads `cifText` into a validated structure. `fileBaseName` (the file name without its extension) is the
    /// display-name fallback when the `data_` block says nothing.
    public static func read(_ cifText: String, fileBaseName: String = "") throws(CIFImportError) -> CIFDocument {
        let parsed = try parse(cifText)
        let cell = UnitCell(a: parsed.a, b: parsed.b, c: parsed.c,
                            alphaDeg: parsed.alpha, betaDeg: parsed.beta, gammaDeg: parsed.gamma)
        let name = displayName(dataBlockName: parsed.dataBlockName, fileBaseName: fileBaseName)
        // The symmetry is detected from the (expanded) cell either way (registration 3.4, amended 2026-10-14): a
        // file's loop may declare less than its atoms have.
        let structure = CrystalStructure.detectingSymmetry(
            cell: cell, sites: parsed.sites, fileOperations: parsed.operations,
            spaceGroupNumber: parsed.spaceGroupNumber, name: name, coarsestHalfStep: parsed.coarsestHalfStep)
        let issues = structure.validationIssues
        guard issues.isEmpty else { throw .invalidStructure(issues.map(\.message)) }
        return CIFDocument(structure: structure, dataBlockName: parsed.dataBlockName, displayName: name,
                           declaredSpaceGroupSymbol: parsed.spaceGroupSymbol,
                           declaredSpaceGroupNumber: parsed.spaceGroupNumber,
                           coarsestHalfStep: parsed.coarsestHalfStep,
                           hadSymmetryOperations: parsed.operations != nil)
    }

    /// The name a structure imports under: the CIF's `data_` block name, unless it says nothing — empty, or the
    /// placeholder a structure editor writes into every file (VESTA exports `data_VESTA_phase_1` whatever the
    /// structure). The file name is then the reader's only handle.
    public static func displayName(dataBlockName: String?, fileBaseName: String) -> String {
        guard let block = dataBlockName?.trimmingCharacters(in: .whitespaces), !block.isEmpty,
              block.range(of: #"^VESTA_phase_[0-9]+$"#, options: .regularExpression) == nil
        else { return fileBaseName }
        return block
    }

    // MARK: - Parsed intermediate structure

    private struct ParsedStructure {
        var dataBlockName: String?
        var a: Double, b: Double, c: Double
        var alpha: Double, beta: Double, gamma: Double
        var sites: [AtomSite]
        /// The expanding operations (deduplicated), or nil when the file has none and the symmetry is to be detected.
        var operations: [SymmetryOperation]?
        var coarsestHalfStep: Double
        var spaceGroupNumber: Int?
        var spaceGroupSymbol: String?
    }

    /// One asymmetric-unit site plus how precisely each of its coordinates was written: `halfStep` is the half-width
    /// of the rounding interval of the literal text, 5e-4 for `0.333`, 0 for `0.2500` (exactly ¼).
    private struct BaseSite {
        let z: Int
        let label: String
        let occupancy: Double
        let position: SIMD3<Double>
        let halfStep: SIMD3<Double>
    }

    private static func parse(_ text: String) throws(CIFImportError) -> ParsedStructure {
        let tokens = tokenize(text)

        var dataBlockName: String?
        var blocks: [(name: String, definesStructure: Bool, hasSymmetryLoop: Bool)] = []
        var cellValues: [String: String] = [:]
        var atomSiteRows: [[String: String]] = []
        var symmetryOps: [String] = []
        var haveSymmetryLoop = false

        func isTag(_ token: String) -> Bool { token.hasPrefix("_") }
        func isLoop(_ token: String) -> Bool { token.lowercased() == "loop_" }
        func isDataBlock(_ token: String) -> Bool { token.lowercased().hasPrefix("data_") }

        var index = 0
        while index < tokens.count {
            let token = tokens[index]

            if isDataBlock(token) {
                let name = String(token.dropFirst("data_".count))
                blocks.append((name: name.isEmpty ? "(unnamed)" : name, definesStructure: false, hasSymmetryLoop: false))
                if dataBlockName == nil {
                    dataBlockName = name.isEmpty ? nil : name
                }
                index += 1
                continue
            }

            if isLoop(token) {
                index += 1
                var loopTags: [String] = []
                while index < tokens.count, isTag(tokens[index]) {
                    loopTags.append(tokens[index].lowercased())
                    index += 1
                }
                guard !loopTags.isEmpty else { continue }

                var values: [String] = []
                while index < tokens.count,
                      !isTag(tokens[index]), !isLoop(tokens[index]), !isDataBlock(tokens[index]) {
                    values.append(tokens[index])
                    index += 1
                }
                // A well-formed loop is always a whole number of rows; a remainder means the file was cut mid-row.
                guard values.count % loopTags.count == 0 else {
                    throw .truncatedLoop(firstTag: loopTags[0], columns: loopTags.count, values: values.count)
                }
                let rowCount = values.count / loopTags.count
                var rows: [[String: String]] = []
                rows.reserveCapacity(rowCount)
                for row in 0..<rowCount {
                    var fields: [String: String] = [:]
                    for (column, tag) in loopTags.enumerated() {
                        fields[tag] = values[row * loopTags.count + column]
                    }
                    rows.append(fields)
                }

                if loopTags.contains("_atom_site_fract_x") {
                    if !blocks.isEmpty { blocks[blocks.count - 1].definesStructure = true }
                    atomSiteRows.append(contentsOf: rows)
                } else if loopTags.contains("_symmetry_equiv_pos_as_xyz")
                    || loopTags.contains("_space_group_symop_operation_xyz") {
                    haveSymmetryLoop = true
                    if !blocks.isEmpty { blocks[blocks.count - 1].hasSymmetryLoop = true }
                    let key = loopTags.contains("_symmetry_equiv_pos_as_xyz")
                        ? "_symmetry_equiv_pos_as_xyz" : "_space_group_symop_operation_xyz"
                    for row in rows {
                        if let op = row[key] { symmetryOps.append(op) }
                    }
                }
                continue
            }

            if isTag(token) {
                let tag = token.lowercased()
                guard index + 1 < tokens.count else { break }
                if tag == "_cell_length_a", !blocks.isEmpty {
                    blocks[blocks.count - 1].definesStructure = true
                }
                cellValues[tag] = tokens[index + 1]
                index += 2
                continue
            }

            // A stray value with no owning tag or loop: skip it and keep scanning.
            index += 1
        }

        let structureBlocks = blocks.filter(\.definesStructure)
        if structureBlocks.count > 1 {
            throw .multipleStructures(blockNames: structureBlocks.map(\.name))
        }
        if !structureBlocks.isEmpty {
            let symmetryOnly = blocks.filter { !$0.definesStructure && $0.hasSymmetryLoop }
            if !symmetryOnly.isEmpty {
                throw .symmetryOperationsOutsideStructure(blockNames: symmetryOnly.map(\.name))
            }
        }

        func requiredNumber(_ tag: String) throws(CIFImportError) -> Double {
            guard let raw = cellValues[tag] else { throw .missingTag(tag) }
            return try parseNumber(raw, tag: tag)
        }

        let a = try requiredNumber("_cell_length_a")
        let b = try requiredNumber("_cell_length_b")
        let c = try requiredNumber("_cell_length_c")
        let alpha = try requiredNumber("_cell_angle_alpha")
        let beta = try requiredNumber("_cell_angle_beta")
        let gamma = try requiredNumber("_cell_angle_gamma")

        var spaceGroupNumber: Int?
        for tag in ["_symmetry_int_tables_number", "_space_group_it_number"] {
            guard let raw = cellValues[tag]?.trimmingCharacters(in: .whitespaces),
                  !raw.isEmpty, raw != "?", raw != "."
            else { continue }
            var text = raw
            if let paren = text.firstIndex(of: "(") { text = String(text[text.startIndex..<paren]) }
            if let value = Int(text) { spaceGroupNumber = value; break }
        }

        guard !atomSiteRows.isEmpty else { throw .noAtomSites }

        var baseSites: [BaseSite] = []
        baseSites.reserveCapacity(atomSiteRows.count)
        for row in atomSiteRows {
            let label = row["_atom_site_label"] ?? "?"
            guard let xRaw = row["_atom_site_fract_x"] else { throw .missingTag("_atom_site_fract_x") }
            guard let yRaw = row["_atom_site_fract_y"] else { throw .missingTag("_atom_site_fract_y") }
            guard let zRaw = row["_atom_site_fract_z"] else { throw .missingTag("_atom_site_fract_z") }
            let x = try parseNumber(xRaw, tag: "_atom_site_fract_x (\(label))")
            let y = try parseNumber(yRaw, tag: "_atom_site_fract_y (\(label))")
            let zCoord = try parseNumber(zRaw, tag: "_atom_site_fract_z (\(label))")
            let occupancy: Double
            if let occRaw = row["_atom_site_occupancy"] {
                occupancy = try parseNumber(occRaw, tag: "_atom_site_occupancy (\(label))")
            } else {
                occupancy = 1.0
            }
            let symbolSource = row["_atom_site_type_symbol"] ?? label
            let symbol = Element.normalizedSymbol(fromCIF: symbolSource)
            guard let atomicNumber = Element.atomicNumber(forSymbol: symbol) else {
                throw .unknownElementSymbol(symbolSource, atomSiteLabel: label)
            }
            baseSites.append(BaseSite(
                z: atomicNumber, label: label, occupancy: occupancy,
                position: SIMD3(x, y, zCoord),
                halfStep: SIMD3(coordinateHalfStep(xRaw), coordinateHalfStep(yRaw), coordinateHalfStep(zRaw))))
        }

        var ops: [SymmetryOperation] = []
        for raw in symmetryOps {
            do { ops.append(try SymmetryOperation(xyz: raw)) } catch { throw .unsupportedSymmetryOperation(raw) }
        }

        let sites: [AtomSite]
        let operations: [SymmetryOperation]?
        if haveSymmetryLoop, ops.contains(where: { !$0.isIdentity }) {
            // A file that names its space group must list every operation of it (in the setting its cell is written
            // in). Fewer builds a smaller structure (rock salt without its F-centring is CsCl, still cubic m-3m).
            if let number = spaceGroupNumber,
               let required = requiredOperationCount(
                   spaceGroup: number, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma) {
                let listed = SymmetryOperation.distinctCount(ops)
                if listed < required {
                    throw .incompleteSymmetryOperations(
                        spaceGroupNumber: number, declaredName: declaredHermannMauguin(cellValues),
                        listed: listed, required: required)
                }
            }
            let lattice = UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma).realLattice
            sites = expand(baseSites, using: ops, latReal: lattice)
            var distinct: [SymmetryOperation] = []
            for op in ops where !distinct.contains(where: { $0.isEquivalent(to: op) }) {
                distinct.append(SymmetryOperation(rotation: op.rotation, translation: op.translationMod1))
            }
            operations = distinct
        } else {
            // No expanding operations. Taking the listed sites as the complete cell is legitimate only when the file
            // itself says they are — P1, or no space-group declaration at all. A file that declares a non-P1 group
            // and omits its operations is an asymmetric unit with no way to expand it. (py4DSTEM hands CIF parsing to
            // pymatgen, which can derive the operations from the symbol; this importer carries no space-group tables,
            // so the honest alternative is refusing.)
            if let declared = declaredNonP1SpaceGroup(cellValues) {
                throw .missingSymmetryOperations(spaceGroup: declared)
            }
            sites = baseSites.map {
                AtomSite(z: $0.z, fractional: wrapFractional($0.position), occupancy: $0.occupancy, label: $0.label)
            }
            operations = nil
        }

        let coarsestHalfStep = baseSites.reduce(0.0) {
            max($0, max($1.halfStep.x, max($1.halfStep.y, $1.halfStep.z)))
        }
        return ParsedStructure(
            dataBlockName: dataBlockName,
            a: a, b: b, c: c, alpha: alpha, beta: beta, gamma: gamma,
            sites: sites, operations: operations, coarsestHalfStep: coarsestHalfStep,
            spaceGroupNumber: spaceGroupNumber, spaceGroupSymbol: declaredHermannMauguin(cellValues))
    }

    // MARK: - Tokenizer

    /// Turns CIF text into a flat token stream: `data_...` block markers, `loop_`, `_tag` names, and (possibly quoted
    /// or semicolon-delimited) values, with `#` comments stripped outside of quotes and text fields.
    private static func tokenize(_ text: String) -> [String] {
        let lines = text.components(separatedBy: .newlines)
        var tokens: [String] = []
        var index = 0
        while index < lines.count {
            let line = lines[index]
            if line.hasPrefix(";") {
                var content = String(line.dropFirst())
                index += 1
                while index < lines.count, !lines[index].hasPrefix(";") {
                    content += "\n" + lines[index]
                    index += 1
                }
                if index < lines.count { index += 1 }
                let trimmed = content.trimmingCharacters(in: .whitespacesAndNewlines)
                if !trimmed.isEmpty { tokens.append(trimmed) }
                continue
            }
            tokens.append(contentsOf: tokenizeLine(line))
            index += 1
        }
        return tokens
    }

    private static func tokenizeLine(_ line: String) -> [String] {
        var result: [String] = []
        let chars = Array(line)
        var i = 0
        while i < chars.count {
            let c = chars[i]
            if c == " " || c == "\t" { i += 1; continue }
            if c == "#" { break }
            if c == "'" || c == "\"" {
                let quote = c
                var j = i + 1
                var value = ""
                while j < chars.count {
                    if chars[j] == quote {
                        let next = j + 1 < chars.count ? chars[j + 1] : " "
                        if next == " " || next == "\t" { break }
                    }
                    value.append(chars[j])
                    j += 1
                }
                result.append(value)
                i = min(j + 1, chars.count)
                continue
            }
            var value = ""
            while i < chars.count, !chars[i].isWhitespace, chars[i] != "#" {
                value.append(chars[i])
                i += 1
            }
            result.append(value)
        }
        return result
    }

    // MARK: - Numbers

    /// Strips a CIF standard-uncertainty suffix, e.g. `4.0782(3)` → `4.0782`.
    private static func parseNumber(_ raw: String, tag: String) throws(CIFImportError) -> Double {
        var text = raw
        if let parenIndex = text.firstIndex(of: "(") {
            text = String(text[text.startIndex..<parenIndex])
        }
        text = text.trimmingCharacters(in: .whitespaces)
        guard let value = Double(text) else { throw .unparsableNumber(tag: tag, value: raw) }
        return value
    }

    /// How far a written fractional coordinate may sit from its true value, judged from the literal text alone:
    /// `10⁻ᵈ` for d decimals (CIF writers truncate as often as they round, so a full ulp of the last decimal), 0 for
    /// a value exactly on a small dyadic rational (`0.2500` is ¼), capped at 1e-2.
    public static func coordinateHalfStep(_ raw: String) -> Double {
        var text = raw
        if let parenIndex = text.firstIndex(of: "(") {
            text = String(text[text.startIndex..<parenIndex])
        }
        text = text.trimmingCharacters(in: .whitespaces)
        guard let value = Double(text) else { return 0 }

        for denominator in [1.0, 2, 4, 8, 16] {
            let scaled = value * denominator
            if abs(scaled - scaled.rounded()) <= 1e-12 * max(1, abs(scaled)) { return 0 }
        }

        guard let dot = text.firstIndex(of: ".") else { return 1e-2 }
        let decimals = text.distance(from: text.index(after: dot), to: text.endIndex)
        return min(1e-2, pow(10, -Double(decimals)))
    }

    // MARK: - Space-group order

    /// Order of the point group of each range of IT numbers (ITA crystal classes).
    private static let pointGroupOrders: [(numbers: ClosedRange<Int>, order: Int)] = [
        (1...1, 1), (2...9, 2), (10...46, 4), (47...74, 8), (75...82, 4),
        (83...122, 8), (123...142, 16), (143...146, 3), (147...161, 6),
        (162...167, 12), (168...174, 6), (175...190, 12), (191...194, 24),
        (195...199, 12), (200...220, 24), (221...230, 48),
    ]

    /// IT numbers of the centred groups, by lattice points per conventional cell: C/A 2, I 2, R (hexagonal axes) 3, F 4.
    private static let centredGroups: [(numbers: Set<Int>, multiplicity: Int)] = [
        (Set([5, 8, 9, 12, 15, 20, 21, 35, 36, 37, 38, 39, 40, 41, 63, 64, 65, 66, 67, 68,
              23, 24, 44, 45, 46, 71, 72, 73, 74, 79, 80, 82, 87, 88, 97, 98,
              107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142,
              197, 199, 204, 206, 211, 214, 217, 220, 229, 230]), 2),
        (Set([146, 148, 155, 160, 161, 166, 167]), 3),
        (Set([22, 42, 43, 69, 70, 196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]), 4),
    ]

    /// Number of general-position operations of space group `number` (1 ... 230) in the conventional setting, centring
    /// translations included: point-group order × lattice points per cell (Fm-3m 48 × 4 = 192, P6₃/mmc 24, R-3m in
    /// hexagonal axes 36). nil outside 1 ... 230.
    public static func spaceGroupOrder(_ number: Int) -> Int? {
        guard let point = pointGroupOrders.first(where: { $0.numbers.contains(number) })?.order
        else { return nil }
        let centring = centredGroups.first(where: { $0.numbers.contains(number) })?.multiplicity ?? 1
        return point * centring
    }

    /// The fewest distinct operations a CIF for `number` may list before it is called incomplete: the conventional
    /// order, except for a centred group written in a cell that is not its conventional one (a primitive cell lists
    /// order / centring; rhombohedral axes a third), judged from the angles alone.
    public static func requiredOperationCount(
        spaceGroup number: Int, alphaDeg: Double, betaDeg: Double, gammaDeg: Double
    ) -> Int? {
        guard let order = spaceGroupOrder(number) else { return nil }
        let centring = centredGroups.first(where: { $0.numbers.contains(number) })?.multiplicity ?? 1
        guard centring > 1 else { return order }
        func right(_ angle: Double) -> Bool { abs(angle - 90) <= 0.05 }
        let conventional: Bool
        switch number {
        case 3...15:
            conventional = [alphaDeg, betaDeg, gammaDeg].filter(right).count >= 2
        case 143...194:
            conventional = right(alphaDeg) && right(betaDeg) && abs(gammaDeg - 120) <= 0.05
        default:
            conventional = right(alphaDeg) && right(betaDeg) && right(gammaDeg)
        }
        return conventional ? order : order / centring
    }

    /// The Hermann–Mauguin symbol the file states; nil when it states none.
    private static func declaredHermannMauguin(_ cellValues: [String: String]) -> String? {
        for tag in ["_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt",
                    "_space_group_name_h-m_full", "_space_group_name_h-m_ref"] {
            if let raw = cellValues[tag]?.trimmingCharacters(in: .whitespaces),
               !raw.isEmpty, raw != "?", raw != "." { return raw }
        }
        return nil
    }

    /// The space group the file itself declares, when that declaration is anything other than P1: a display name for
    /// the refusal, or nil when the file declares P1 or declares nothing (`?` and `.` count as absent). A declaration
    /// that exists but cannot be read is returned too: guessing P1 from an unreadable tag would be the silent
    /// assumption this guard removes.
    private static func declaredNonP1SpaceGroup(_ cellValues: [String: String]) -> String? {
        func nonNull(_ tag: String) -> String? {
            guard let raw = cellValues[tag]?.trimmingCharacters(in: .whitespaces),
                  !raw.isEmpty, raw != "?", raw != "." else { return nil }
            return raw
        }
        for tag in ["_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt",
                    "_space_group_name_h-m_full", "_space_group_name_h-m_ref",
                    "_space_group_name_hall", "_symmetry_space_group_name_hall"] {
            guard let raw = nonNull(tag) else { continue }
            let normalized = raw.lowercased().filter { !$0.isWhitespace }
            if normalized != "p1" { return raw }
        }
        for tag in ["_space_group_it_number", "_symmetry_int_tables_number"] {
            guard let raw = nonNull(tag) else { continue }
            guard let number = Int(raw) else { return "\(tag) = \(raw)" }
            if number != 1 { return "No. \(number)" }
        }
        return nil
    }

    // MARK: - Expansion

    /// Applies every operation to every asymmetric-unit site, wraps into [0, 1), and merges images that land on the
    /// same site: within the rounding error the written precision allows for both of them (propagated through the
    /// operation's rows), and within `maxMergeDistance` in real space, so a coarse cell cannot collapse chemistry. A
    /// pair too far apart to be rounding but too close to be chemistry is left in place and refused by the
    /// close-contact rule rather than silently deleted.
    private static func expand(
        _ baseSites: [BaseSite], using ops: [SymmetryOperation], latReal: [SIMD3<Double>]
    ) -> [AtomSite] {
        let maxMergeDistance = 0.10
        let usableLattice = latReal.allSatisfy { $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }

        struct Image {
            let row: Int
            let z: Int
            let label: String
            let occupancy: Double
            let position: SIMD3<Double>
            let uncertainty: SIMD3<Double>
        }

        var kept: [Image] = []
        for (row, site) in baseSites.enumerated() {
            for op in ops {
                let mapped = wrapFractional(op.apply(site.position))
                let uncertainty = SIMD3<Double>(
                    propagate(op.rotation.rows[0], site.halfStep),
                    propagate(op.rotation.rows[1], site.halfStep),
                    propagate(op.rotation.rows[2], site.halfStep))
                // Only images of the SAME row merge (Gate B amendment, 2026-10-15): two rows on one position are two
                // species sharing a site (Fe2+ and Fe3+ in magnetite), each with its own occupancy; at full
                // occupancy they are refused by the close-contact rule, never silently merged.
                let isDuplicate = kept.contains { existing in
                    guard existing.row == row else { return false }
                    let tolerance = existing.uncertainty + uncertainty
                    guard fractionalAxisDelta(existing.position.x, mapped.x) <= max(tolerance.x, 1e-6),
                          fractionalAxisDelta(existing.position.y, mapped.y) <= max(tolerance.y, 1e-6),
                          fractionalAxisDelta(existing.position.z, mapped.z) <= max(tolerance.z, 1e-6)
                    else { return false }
                    guard usableLattice else { return true }
                    return CrystalStructure.minimumImageDistanceSquared(existing.position, mapped, latReal)
                        .squareRoot() <= maxMergeDistance
                }
                if !isDuplicate {
                    kept.append(Image(row: row, z: site.z, label: site.label, occupancy: site.occupancy,
                                      position: mapped, uncertainty: uncertainty))
                }
            }
        }
        return kept.map { AtomSite(z: $0.z, fractional: $0.position, occupancy: $0.occupancy, label: $0.label) }
    }

    /// Rounding error of one output coordinate: `|c_x|·δx + |c_y|·δy + |c_z|·δz` for the operator row.
    private static func propagate(_ row: SIMD3<Int>, _ halfStep: SIMD3<Double>) -> Double {
        Double(abs(row.x)) * halfStep.x + Double(abs(row.y)) * halfStep.y + Double(abs(row.z)) * halfStep.z
    }
}
