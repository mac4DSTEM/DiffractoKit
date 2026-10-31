//
//  CIFImport.swift
//  Role: the subset of CIF that defines a structure (cell, the atom-site loop, the symmetry-operation loop), read into
//        a `CrystalStructure`. Ported from reference commit 7448401b `CIFImport.swift`: the
//        tokenizer, the `(3)` uncertainty suffix, the element symbol rule, the written-precision half-step, the
//        expansion that merges images within the propagated half-step AND 0.10 Å, and the refusals that protect the
//        structure. Generalised: every crystal system is admitted; a file without operations has its symmetry
//        detected from the sites; the two app-specific refusals (point group, scattering factor) are gone.
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
    /// As ``incompleteSymmetryOperations(spaceGroupNumber:declaredName:listed:required:)`` for a file that names its group
    /// by Hermann–Mauguin symbol alone (no IT number): the order is read from the symbol (``requiredOperationCount(hermannMauguin:alphaDeg:betaDeg:gammaDeg:rightAngleToleranceDeg:rotationsDiagonal:)``).
    /// A file that names its group by Hall symbol alone (no IT number, no readable Hermann–Mauguin symbol) is held to the
    /// order of the group the Hall symbol states (``HallSymbol/operations(_:)``), and `symbol` is the Hall symbol.
    case incompleteSymmetryOperationsForSymbol(symbol: String, listed: Int, required: Int)
    /// A `data_` block that defines no structure carries a symmetry-operation loop while another block defines one.
    case symmetryOperationsOutsideStructure(blockNames: [String])
    /// The operations are enough in number but cannot be the file's space group: one does not map the cell onto itself,
    /// the list is not closed, or its crystal class is not the IT number's (``SymmetryListProblem``).
    case inconsistentSymmetryOperations(SymmetryListProblem)

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
        case .incompleteSymmetryOperationsForSymbol(let symbol, let listed, let required):
            return "The CIF lists \(listed) distinct symmetry operations but space group \"\(symbol)\" has \(required) in its "
                + "conventional setting — the list is incomplete, so expanding the atom sites with it would build a "
                + "smaller structure than the file describes. Re-export the CIF with its full symmetry-operation loop, "
                + "or supply a P1 file that lists every atom in the cell."
        case .symmetryOperationsOutsideStructure(let names):
            return "The CIF has a symmetry-operation loop in a data_ block that holds no structure "
                + "(\(names.joined(separator: ", "))). Reading the file as one crystal would apply those operations to "
                + "another block's atoms. Move the loop into the structure's own data_ block, or split the file."
        case .inconsistentSymmetryOperations(let problem):
            return "The CIF's symmetry operations cannot be its space group: \(problem). Expanding the atom sites with "
                + "them would build a wrong structure, so the file is refused. Re-export the CIF with the operation "
                + "loop of its own cell and setting, or supply a P1 file that lists every atom in the cell."
        }
    }
}

/// `localizedDescription` (what an app shows) is the same text as `description`.
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
    /// The Hall symbol the file states, as written; nil when none.
    public let declaredHallSymbol: String?
    /// The coarsest rounding half-step over every written fractional coordinate (0 when all are exact).
    public let coarsestHalfStep: Double
    /// Whether the file carried expanding operations (`structure.fileOperations`): a symmetry-operation loop, or, with no
    /// loop, a Hall symbol read to its group. The structure's operations are detected from the expanded cell either way.
    public let hadSymmetryOperations: Bool
}

public enum CIFImport {

    /// The thresholds of the import that are the suite's choices rather than the CIF grammar (DiffractoKit 0.2.0), with
    /// `Policy.default` the set every release before used. A distance, cap or bound that is NaN or negative is taken as
    /// 0 (merge nothing but coincident images, treat coordinates as exact, parse no number but 0); so are the two metric
    /// tolerances, and an infinite one too (+infinity is read as 0, not as "accept everything"). The position tolerance
    /// is not clamped: one outside its domain (`SymmetryDetection.Tolerance`) makes `read` throw `invalidStructure`; the
    /// tolerances of the nested policies are described where they are declared.
    public struct Policy: Sendable, Equatable {
        /// Images of one site merge only within this distance (Å) as well as within the propagated rounding error, so a
        /// coarse cell cannot collapse chemistry. Default 0.10.
        public var maximumMergeDistance: Double
        /// The floor under the propagated rounding error of the merge (fractional). Default 1e-6.
        public var minimumMergeTolerance: Double
        /// The cap on a written coordinate's rounding half-step (fractional). Default 1e-2.
        public var maximumCoordinateHalfStep: Double
        /// A cell angle within this of 90 degrees (or 120 for hexagonal axes) is called right when the centring of a
        /// declared group is judged from the angles. Default 0.05.
        public var rightAngleToleranceDeg: Double
        /// The largest number a symmetry operation's text may hold. Default 1000.
        public var maximumOperationMagnitude: Double
        /// How far a listed rotation may change the cell metric and still count as a symmetry of the cell: each entry of
        /// |RᵀGR − G| over sqrt(G_ii G_jj). Default 0.01 (about 0.6 degrees in an angle, 0.5 % in a length; an operation
        /// of another cell or setting changes it by tens of percent). Only this check of the operation list has a knob:
        /// closure and the crystal class are not tolerances and cannot be turned off.
        public var operationMetricTolerance: Double
        /// How two positions are called the same site in the symmetry detection. Default fractional 5e-3.
        public var positionTolerance: SymmetryDetection.Tolerance
        /// The metric family tolerances (relative on lengths, degrees on angles). Default 1e-9 both.
        public var metricLengthTolerance: Double
        public var metricAngleToleranceDeg: Double
        /// The occupancy tolerance, translation grid and precision-limit factor of the symmetry detection.
        public var symmetry: SymmetryDetection.Policy
        /// The close-contact and split-site thresholds of the validation the import ends with.
        public var validation: CrystalStructure.ValidationPolicy

        public init(maximumMergeDistance: Double = 0.10, minimumMergeTolerance: Double = 1e-6,
                    maximumCoordinateHalfStep: Double = 1e-2, rightAngleToleranceDeg: Double = 0.05,
                    maximumOperationMagnitude: Double = SymmetryOperation.defaultMaximumMagnitude,
                    operationMetricTolerance: Double = 0.01,
                    positionTolerance: SymmetryDetection.Tolerance = .default,
                    metricLengthTolerance: Double = SymmetryDetection.defaultMetricLengthTolerance,
                    metricAngleToleranceDeg: Double = SymmetryDetection.defaultMetricAngleToleranceDeg,
                    symmetry: SymmetryDetection.Policy = .default,
                    validation: CrystalStructure.ValidationPolicy = .default) {
            self.maximumMergeDistance = maximumMergeDistance
            self.minimumMergeTolerance = minimumMergeTolerance
            self.maximumCoordinateHalfStep = maximumCoordinateHalfStep
            self.rightAngleToleranceDeg = rightAngleToleranceDeg
            self.maximumOperationMagnitude = maximumOperationMagnitude
            self.operationMetricTolerance = operationMetricTolerance
            self.positionTolerance = positionTolerance
            self.metricLengthTolerance = metricLengthTolerance
            self.metricAngleToleranceDeg = metricAngleToleranceDeg
            self.symmetry = symmetry
            self.validation = validation
        }

        /// The suite's values.
        public static let `default` = Policy()

        func nonNegative(_ value: Double) -> Double { value.isFinite && value >= 0 ? value : 0 }
    }

    /// Reads `cifText` into a validated structure. `fileBaseName` (the file name without its extension) is the
    /// display-name fallback when the `data_` block says nothing.
    public static func read(_ cifText: String, fileBaseName: String = "",
                            policy: Policy = .default) throws(CIFImportError) -> CIFDocument {
        let parsed = try parse(cifText, policy: policy)
        let cell = UnitCell(a: parsed.a, b: parsed.b, c: parsed.c,
                            alphaDeg: parsed.alpha, betaDeg: parsed.beta, gammaDeg: parsed.gamma)
        let name = displayName(dataBlockName: parsed.dataBlockName, fileBaseName: fileBaseName)
        // The symmetry is detected from the (expanded) cell either way: a
        // file's loop may declare less than its atoms have.
        let structure = CrystalStructure.detectingSymmetry(
            cell: cell, sites: parsed.sites, fileOperations: parsed.operations,
            spaceGroupNumber: parsed.spaceGroupNumber, name: name, coarsestHalfStep: parsed.coarsestHalfStep,
            metricLengthTolerance: policy.nonNegative(policy.metricLengthTolerance),
            metricAngleToleranceDeg: policy.nonNegative(policy.metricAngleToleranceDeg),
            tolerance: policy.positionTolerance, policy: policy.symmetry)
        let issues = structure.validationIssues(using: policy.validation)
        guard issues.isEmpty else { throw .invalidStructure(issues.map(\.message)) }
        return CIFDocument(structure: structure, dataBlockName: parsed.dataBlockName, displayName: name,
                           declaredSpaceGroupSymbol: parsed.spaceGroupSymbol,
                           declaredSpaceGroupNumber: parsed.spaceGroupNumber,
                           declaredHallSymbol: parsed.hallSymbol,
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
        var hallSymbol: String?
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

    private static func parse(_ text: String, policy: Policy) throws(CIFImportError) -> ParsedStructure {
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
                halfStep: SIMD3(coordinateHalfStep(xRaw, cap: policy.maximumCoordinateHalfStep),
                                coordinateHalfStep(yRaw, cap: policy.maximumCoordinateHalfStep),
                                coordinateHalfStep(zRaw, cap: policy.maximumCoordinateHalfStep))))
        }

        var ops: [SymmetryOperation] = []
        for raw in symmetryOps {
            do { ops.append(try SymmetryOperation(xyz: raw, maximumMagnitude: policy.maximumOperationMagnitude)) } catch { throw .unsupportedSymmetryOperation(raw) }
        }

        // A Hall symbol states the whole group in the file's own setting, so a file that names its group by Hall symbol
        // and carries no loop is expanded with the symbol's operations (`HallSymbol`, exact against spglib's 530
        // settings). pymatgen means its symbol fallback to do this, but looks a Hall symbol up among Hermann–Mauguin
        // names, finds none and expands as P1 with a warning (pymatgen-core 2026.10.16, `CifParser.get_symops`).
        let fileListsOperations = haveSymmetryLoop && ops.contains(where: { !$0.isIdentity })
        let hallSymbol = declaredHall(cellValues)
        let hallOperations = hallSymbol.flatMap { symbol -> [SymmetryOperation]? in try? HallSymbol.operations(symbol) }
        // A loop, even one holding only the identity, is what the file lists and wins over the symbol (as in pymatgen).
        let expandsWithHall = !haveSymmetryLoop && (hallOperations?.count ?? 0) > 1
        if expandsWithHall, let hallOperations { ops = hallOperations }

        let sites: [AtomSite]
        let operations: [SymmetryOperation]?
        if fileListsOperations || expandsWithHall {
            // A file that names its space group must list every operation of it (in the setting its cell is written
            // in). Fewer builds a smaller structure (rock salt without its F-centring is CsCl, still cubic m-3m).
            // The judgement of a monoclinic centred cell (conventional or primitive) also reads the listed rotations:
            // angles alone cannot tell a primitive C-cell from a conventional one when beta = 90.
            let diagonal = ops.allSatisfy { op in
                (0..<3).allSatisfy { r in (0..<3).allSatisfy { c in r == c || op.rotation[r, c] == 0 } }
            }
            let listed = SymmetryOperation.distinctCount(ops)
            // A Hall group stands in for the loop here: a Hall symbol smaller than the IT number or Hermann–Mauguin
            // symbol the file also states (a mistyped "-P 1" beside No. 225) is refused like a short loop.
            if let number = spaceGroupNumber,
               let required = requiredOperationCount(
                   spaceGroup: number, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma,
                   rightAngleToleranceDeg: policy.rightAngleToleranceDeg, rotationsDiagonal: diagonal) {
                if listed < required {
                    throw .incompleteSymmetryOperations(
                        spaceGroupNumber: number, declaredName: declaredHermannMauguin(cellValues),
                        listed: listed, required: required)
                }
            } else if spaceGroupNumber == nil, let symbol = declaredHermannMauguin(cellValues),
                      let required = requiredOperationCount(
                          hermannMauguin: symbol, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma,
                          rightAngleToleranceDeg: policy.rightAngleToleranceDeg, rotationsDiagonal: diagonal),
                      listed < required {
                throw .incompleteSymmetryOperationsForSymbol(symbol: symbol, listed: listed, required: required)
            } else if fileListsOperations, spaceGroupNumber == nil,
                      declaredHermannMauguin(cellValues).flatMap(hermannMauguinFacts) == nil,
                      let hallSymbol, let hallOperations, listed < hallOperations.count {
                // Named by Hall symbol alone: the symbol's group, centring included, is the order the loop must reach.
                throw .incompleteSymmetryOperationsForSymbol(symbol: hallSymbol, listed: listed,
                                                             required: hallOperations.count)
            }
            // A count cannot see wrong operations: each must map the cell onto itself, the list must be closed, and
            // its crystal class must be the IT number's (a Hall group is closed, but its setting may not fit the cell).
            let unitCell = UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma)
            if let problem = operationListProblem(ops, cell: unitCell, spaceGroupNumber: spaceGroupNumber,
                                                  metricTolerance: policy.operationMetricTolerance) {
                throw .inconsistentSymmetryOperations(problem)
            }
            let lattice = unitCell.realLattice
            sites = expand(baseSites, using: ops, latReal: lattice, policy: policy)
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
            // so the honest alternative is refusing. A readable Hall symbol needs no table and was used above.)
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
            spaceGroupNumber: spaceGroupNumber, spaceGroupSymbol: declaredHermannMauguin(cellValues),
            hallSymbol: hallSymbol)
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
    public static func coordinateHalfStep(_ raw: String, cap: Double = 1e-2) -> Double {
        let ceiling: Double = cap.isFinite && cap >= 0 ? cap : 0
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

        guard let dot = text.firstIndex(of: ".") else { return ceiling }
        let decimals = text.distance(from: text.index(after: dot), to: text.endIndex)
        return min(ceiling, pow(10, -Double(decimals)))
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
    ///
    /// `rotationsDiagonal` (the listed rotations all diagonal in the file's basis) refines the monoclinic judgement: a
    /// primitive C-centred monoclinic cell with beta = 90 has two right angles like the conventional one, but its
    /// rotations are not diagonal. nil judges from the angles alone, as before.
    public static func requiredOperationCount(
        spaceGroup number: Int, alphaDeg: Double, betaDeg: Double, gammaDeg: Double,
        rightAngleToleranceDeg: Double = 0.05, rotationsDiagonal: Bool? = nil
    ) -> Int? {
        guard let order = spaceGroupOrder(number) else { return nil }
        let centring = centredGroups.first(where: { $0.numbers.contains(number) })?.multiplicity ?? 1
        guard centring > 1 else { return order }
        let angleBound: Double = rightAngleToleranceDeg.isFinite && rightAngleToleranceDeg >= 0 ? rightAngleToleranceDeg : 0
        func right(_ angle: Double) -> Bool { abs(angle - 90) <= angleBound }
        let conventional: Bool
        switch number {
        case 3...15:
            conventional = [alphaDeg, betaDeg, gammaDeg].filter(right).count >= 2 && rotationsDiagonal != false
        case 143...194:
            conventional = right(alphaDeg) && right(betaDeg) && abs(gammaDeg - 120) <= angleBound
        default:
            conventional = right(alphaDeg) && right(betaDeg) && right(gammaDeg)
        }
        return conventional ? order : order / centring
    }

    /// The point-group orders of the 32 crystal classes by short symbol, with the crystal system's first letter
    /// (a anorthic, m monoclinic, o orthorhombic, t tetragonal, r trigonal, h hexagonal, c cubic).
    private static let pointGroupBySymbol: [String: (order: Int, system: Character)] = [
        "1": (1, "a"), "-1": (2, "a"),
        "2": (2, "m"), "m": (2, "m"), "2/m": (4, "m"),
        "222": (4, "o"), "mm2": (4, "o"), "mmm": (8, "o"),
        "4": (4, "t"), "-4": (4, "t"), "4/m": (8, "t"), "422": (8, "t"), "4mm": (8, "t"), "-42m": (8, "t"),
        "-4m2": (8, "t"), "4/mmm": (16, "t"),
        "3": (3, "r"), "-3": (6, "r"), "32": (6, "r"), "3m": (6, "r"), "-3m": (12, "r"),
        "6": (6, "h"), "-6": (6, "h"), "6/m": (12, "h"), "622": (12, "h"), "6mm": (12, "h"), "-6m2": (12, "h"),
        "-62m": (12, "h"), "6/mmm": (24, "h"),
        "23": (12, "c"), "m-3": (24, "c"), "432": (24, "c"), "-43m": (24, "c"), "m-3m": (48, "c"),
    ]

    /// The lattice letter, the point-group order and the crystal system a Hermann–Mauguin symbol states (screw axes read as
    /// their rotation, glides as mirrors, the placeholder 1 of a monoclinic or trigonal symbol dropped), or nil when the
    /// symbol is not one of the 32 classes (a Hall symbol, a typo). A bare "32" token reads as the screw axis 3₂, which
    /// can only under-count the group (a file is then not refused), never over-count it.
    package static func hermannMauguinFacts(_ symbol: String) -> (centring: Int, order: Int, system: Character)? {
        let text = symbol.trimmingCharacters(in: .whitespaces)
        guard let lattice = text.first, "PABCIFRH".contains(lattice) else { return nil }
        let rest = String(text.dropFirst())
        let screws: Set<String> = ["21", "31", "32", "41", "42", "43", "61", "62", "63", "64", "65"]
        let bareScrews = !rest.contains("_") && lattice != "R" && lattice != "H"
        var axes: [(axis: String, mirror: Bool)] = []
        for chunk in rest.split(whereSeparator: { $0.isWhitespace }) {
            let t = Array(chunk)
            var i = 0
            while i < t.count {
                let c = t[i]
                if c == "-" || c.isNumber {
                    var axis = ""
                    if c == "-" { axis = "-"; i += 1 }
                    guard i < t.count, t[i].isNumber else { return nil }
                    axis.append(t[i]); i += 1
                    if i < t.count, t[i] == "_" {
                        i += 1
                        guard i < t.count, t[i].isNumber else { return nil }
                        i += 1
                    } else if i < t.count, t[i].isNumber, bareScrews, i == 1, axis.count == 1, i + 1 == t.count || t[i + 1] == "/",
                              screws.contains(axis + String(t[i])) {
                        i += 1
                    }
                    var mirror = false
                    if i < t.count, t[i] == "/" {
                        i += 1
                        guard i < t.count, t[i].isLetter else { return nil }
                        mirror = true; i += 1
                    }
                    if axis != "1" || mirror { axes.append((axis, mirror)) }
                } else if "mabcndeg".contains(c) {
                    axes.append(("m", false)); i += 1
                } else {
                    return nil
                }
            }
        }
        // The full symbol writes a mirror as "2/m", "4/m": past the first position, or before a cubic 3 or -3, the axis
        // is the short symbol's m (4/m 2/m 2/m is 4/mmm, 2/m -3 is m-3).
        let cubic = axes.contains { $0.axis == "3" || $0.axis == "-3" } && axes.count > 1
        var point = ""
        for (k, a) in axes.enumerated() {
            if a.mirror, axes.count > 1, k > 0 || cubic || a.axis == "2" { point += "m" }
            else { point += a.axis + (a.mirror ? "/m" : "") }
        }
        if point.isEmpty { point = "1" }
        let centring = lattice == "P" ? 1 : lattice == "F" ? 4 : (lattice == "R" || lattice == "H") ? 3 : 2
        if let hit = pointGroupBySymbol[point] { return (centring, hit.order, hit.system) }
        // Orthorhombic symbols permute (P n m a is mmm, P m 2 m is mm2).
        let key = String(point.sorted())
        for (name, hit) in pointGroupBySymbol where hit.system == "o" && String(name.sorted()) == key {
            return (centring, hit.order, hit.system)
        }
        return nil
    }

    /// As ``requiredOperationCount(spaceGroup:alphaDeg:betaDeg:gammaDeg:rightAngleToleranceDeg:rotationsDiagonal:)`` for a
    /// file that names its group by Hermann–Mauguin symbol alone: the point-group order times the lattice points per cell
    /// the symbol states, divided by the centring when the cell is not the conventional one. nil when the symbol is not
    /// one of the 32 classes' (such a file stays unchecked, as before).
    public static func requiredOperationCount(
        hermannMauguin symbol: String, alphaDeg: Double, betaDeg: Double, gammaDeg: Double,
        rightAngleToleranceDeg: Double = 0.05, rotationsDiagonal: Bool? = nil
    ) -> Int? {
        guard let facts = hermannMauguinFacts(symbol) else { return nil }
        let order = facts.order * facts.centring
        guard facts.centring > 1 else { return order }
        let angleBound: Double = rightAngleToleranceDeg.isFinite && rightAngleToleranceDeg >= 0 ? rightAngleToleranceDeg : 0
        func right(_ angle: Double) -> Bool { abs(angle - 90) <= angleBound }
        let conventional: Bool
        switch facts.system {
        case "m": conventional = [alphaDeg, betaDeg, gammaDeg].filter(right).count >= 2 && rotationsDiagonal != false
        case "r", "h": conventional = right(alphaDeg) && right(betaDeg) && abs(gammaDeg - 120) <= angleBound
        default: conventional = right(alphaDeg) && right(betaDeg) && right(gammaDeg)
        }
        return conventional ? order : order / facts.centring
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

    /// The Hall symbol the file states (`_space_group_name_Hall`, or the older `_symmetry_space_group_name_Hall`); nil when
    /// none.
    private static func declaredHall(_ cellValues: [String: String]) -> String? {
        for tag in ["_space_group_name_hall", "_symmetry_space_group_name_hall"] {
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
        _ baseSites: [BaseSite], using ops: [SymmetryOperation], latReal: [SIMD3<Double>], policy: Policy
    ) -> [AtomSite] {
        let maxMergeDistance = policy.nonNegative(policy.maximumMergeDistance)
        let mergeFloor = policy.nonNegative(policy.minimumMergeTolerance)
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
                // Only images of the SAME row merge: two rows on one position are two
                // species sharing a site (Fe2+ and Fe3+ in magnetite), each with its own occupancy; at full
                // occupancy they are refused by the close-contact rule, never silently merged.
                let isDuplicate = kept.contains { existing in
                    guard existing.row == row else { return false }
                    let tolerance = existing.uncertainty + uncertainty
                    guard fractionalAxisDelta(existing.position.x, mapped.x) <= max(tolerance.x, mergeFloor),
                          fractionalAxisDelta(existing.position.y, mapped.y) <= max(tolerance.y, mergeFloor),
                          fractionalAxisDelta(existing.position.z, mapped.z) <= max(tolerance.z, mergeFloor)
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
