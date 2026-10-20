//
//  MaterialsProject.swift
//  Role: the Materials Project (MP) `GET /materials/summary/` response as value types, its decoding, the refusals, the
//        conversion of one document to a `CrystalStructure` (through `CellStandardisation`), and the check of a fetched
//        document against the phase the user expected. Ported from mac4DSTEM `MaterialsProjectImport.swift` (the wire
//        types, `decode`, `firstDocument`, `request`, `check`, the formula parser, and the document-to-structure route
//        of `crystalModel(from:fetchedAt:)`); the model assembly (`PhaseModel`, DiffractoOrientation) and the network
//        fetch stay out. No network access happens here. The wire shape was verified against the MP openapi.json
//        (2026-09-21): `{"data":[SummaryDoc...],"errors":[{"code":int,"message":str}],"meta":{...}}`; unknown keys are
//        ignored, as for CIF tags the reader does not model.
//
//  Symmetry is DETECTED from the atoms (`CrystalStructure.detectingSymmetry`); MP's own `symmetry.crystal_system`
//  and `symbol` are used only to standardise a centred primitive cell and are recorded, never trusted to classify.
//  A site with several species becomes one `AtomSite` per species at the same position (the CIF reader's rule).
//

import Foundation
import simd

// MARK: - Wire types

/// `GET /materials/summary/` response envelope.
public struct MaterialsProjectSummaryResponse: Decodable, Sendable {
    public let data: [MaterialsProjectDocument]
    /// Absent, not `[]`, on a successful live response (mp-134, 2026-09-21); the schema says nullable.
    public let errors: [MaterialsProjectError]

    public init(data: [MaterialsProjectDocument], errors: [MaterialsProjectError]) {
        self.data = data
        self.errors = errors
    }

    private enum CodingKeys: String, CodingKey { case data, errors }

    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        data = try c.decodeIfPresent([MaterialsProjectDocument].self, forKey: .data) ?? []
        errors = try c.decodeIfPresent([MaterialsProjectError].self, forKey: .errors) ?? []
    }
}

/// One `SummaryDoc`. Only the fields used are modelled.
public struct MaterialsProjectDocument: Decodable, Sendable {
    public let materialID: String
    public let formulaPretty: String?
    public let deprecated: Bool?
    public let lastUpdated: String?
    public let structure: MaterialsProjectStructure?
    public let symmetry: MaterialsProjectSymmetry?

    private enum CodingKeys: String, CodingKey {
        case materialID = "material_id"
        case formulaPretty = "formula_pretty"
        case deprecated
        case lastUpdated = "last_updated"
        case structure
        case symmetry
    }

    public init(materialID: String, formulaPretty: String? = nil, deprecated: Bool? = nil, lastUpdated: String? = nil,
                structure: MaterialsProjectStructure? = nil, symmetry: MaterialsProjectSymmetry? = nil) {
        self.materialID = materialID
        self.formulaPretty = formulaPretty
        self.deprecated = deprecated
        self.lastUpdated = lastUpdated
        self.structure = structure
        self.symmetry = symmetry
    }
}

public struct MaterialsProjectStructure: Decodable, Sendable {
    public let lattice: MaterialsProjectLattice?
    public let sites: [MaterialsProjectSite]?

    public init(lattice: MaterialsProjectLattice? = nil, sites: [MaterialsProjectSite]? = nil) {
        self.lattice = lattice
        self.sites = sites
    }
}

public struct MaterialsProjectLattice: Decodable, Sendable {
    public let a: Double
    public let b: Double
    public let c: Double
    public let alpha: Double
    public let beta: Double
    public let gamma: Double
    /// The served Cartesian lattice matrix (rows = a, b, c vectors), when present. Used directly by the conversion; rebuilt
    /// from a, b, c, alpha, beta, gamma (as `UnitCell` does) when absent.
    public let matrix: [[Double]]?

    public init(a: Double, b: Double, c: Double, alpha: Double, beta: Double, gamma: Double, matrix: [[Double]]? = nil) {
        self.a = a; self.b = b; self.c = c
        self.alpha = alpha; self.beta = beta; self.gamma = gamma
        self.matrix = matrix
    }
}

public struct MaterialsProjectSite: Decodable, Sendable {
    public let label: String?
    public let species: [MaterialsProjectSpecies]?
    public let abc: [Double]?

    public init(label: String? = nil, species: [MaterialsProjectSpecies]? = nil, abc: [Double]? = nil) {
        self.label = label
        self.species = species
        self.abc = abc
    }
}

public struct MaterialsProjectSpecies: Decodable, Sendable {
    public let element: String
    public let occu: Double

    public init(element: String, occu: Double) {
        self.element = element
        self.occu = occu
    }
}

public struct MaterialsProjectSymmetry: Decodable, Sendable {
    public let crystalSystem: String?
    public let symbol: String?
    public let number: Int?
    public let pointGroup: String?
    public let hall: String?

    private enum CodingKeys: String, CodingKey {
        case crystalSystem = "crystal_system"
        case symbol
        case number
        case pointGroup = "point_group"
        case hall
    }

    public init(crystalSystem: String? = nil, symbol: String? = nil, number: Int? = nil,
                pointGroup: String? = nil, hall: String? = nil) {
        self.crystalSystem = crystalSystem
        self.symbol = symbol
        self.number = number
        self.pointGroup = pointGroup
        self.hall = hall
    }
}

public struct MaterialsProjectError: Decodable, Sendable {
    public let code: Int
    public let message: String

    public init(code: Int, message: String) {
        self.code = code
        self.message = message
    }
}

// MARK: - Phase expectation

/// What the user expected the fetched material to be, before it is trusted. Both fields are independently optional;
/// `MaterialsProjectImport.check` grades whatever the caller knows.
public struct PhaseExpectation: Equatable, Sendable {
    public var formula: String?
    public var spaceGroupNumber: Int?

    public init(formula: String? = nil, spaceGroupNumber: Int? = nil) {
        self.formula = formula
        self.spaceGroupNumber = spaceGroupNumber
    }
}

public enum PhaseExpectationVerdict: Equatable, Sendable {
    case matches
    /// Neither a formula nor a space-group number was supplied: nothing to grade against.
    case unchecked
    case mismatch(reason: String)
}

// MARK: - Importer

public enum MaterialsProjectImport {

    /// Every named refusal. A deprecated MP entry still imports (the flag is the caller's to record).
    public enum Failure: Error, Equatable, LocalizedError, Sendable {
        /// The material id does not match `^mp-[0-9]+$`.
        case invalidMaterialID(String)
        /// The response's `data` array was empty.
        case emptyResult(String)
        /// The response's `errors` array was non-empty.
        case apiError(code: Int, message: String)
        /// The document is missing its structure, lattice or sites, or carries a coordinate or element that cannot be used.
        case malformed(String)
        /// A served centred primitive cell could not be verified as safely convertible to the conventional cell: silently
        /// classifying an unconverted primitive cell (an fcc primitive is a = b = c, 60 degrees, which reads as no
        /// cubic metric) would mis-score its family and index zone axes wrongly, so it is refused.
        case nonStandardCell(String)

        public var errorDescription: String? {
            switch self {
            case .invalidMaterialID(let id):
                return "\"\(id)\" is not a valid Materials Project material id (expected \"mp-<digits>\")."
            case .emptyResult(let reason):
                return "Materials Project returned no usable result: \(reason)"
            case .apiError(let code, let message):
                return "Materials Project API error \(code): \(message)"
            case .malformed(let reason):
                return "Materials Project response could not be used: \(reason)"
            case .nonStandardCell(let reason):
                return "Materials Project cell could not be safely standardised: \(reason)"
            }
        }
    }

    /// A document converted to a structure.
    public struct ImportedStructure: Sendable, Equatable {
        public let structure: CrystalStructure
        /// Whether the served cell had to be converted from a centred primitive cell.
        public let cell: CellStandardisation.Outcome
    }

    /// The `GET /materials/summary/` request for one material id. The id is validated first: a malformed id comes back
    /// from the API as an empty `data` array, indistinguishable from "not found".
    public static func request(materialID: String, apiKey: String) throws -> URLRequest {
        guard isValidMaterialID(materialID) else { throw Failure.invalidMaterialID(materialID) }
        var components = URLComponents(string: "https://api.materialsproject.org/materials/summary/")!
        components.queryItems = [
            URLQueryItem(name: "material_ids", value: materialID),
            URLQueryItem(name: "_fields", value: "material_id,formula_pretty,structure,symmetry,deprecated,last_updated"),
        ]
        guard let url = components.url else { throw Failure.invalidMaterialID(materialID) }
        var request = URLRequest(url: url)
        request.httpMethod = "GET"
        request.setValue(apiKey, forHTTPHeaderField: "X-API-KEY")
        return request
    }

    private static func isValidMaterialID(_ id: String) -> Bool {
        guard id.hasPrefix("mp-") else { return false }
        let digits = id.dropFirst(3)
        return !digits.isEmpty && digits.allSatisfy { $0.isASCII && $0.isNumber }
    }

    /// Decodes the raw response body. A body that is not the documented envelope throws `DecodingError`.
    public static func decode(_ data: Data) throws -> MaterialsProjectSummaryResponse {
        try JSONDecoder().decode(MaterialsProjectSummaryResponse.self, from: data)
    }

    /// The one document to act on. `errors` is checked before `data`: an API error beside an empty result is the error.
    public static func firstDocument(in response: MaterialsProjectSummaryResponse) throws -> MaterialsProjectDocument {
        if let first = response.errors.first { throw Failure.apiError(code: first.code, message: first.message) }
        guard let doc = response.data.first else { throw Failure.emptyResult("the response's data array was empty.") }
        return doc
    }

    /// The structure of one fetched document, converted to the conventional cell when MP served a centred primitive
    /// cell (`CellStandardisation`), with the symmetry detected from the atoms. Structural absence, an unusable
    /// coordinate or an unknown element throws `Failure.malformed`; a conversion that fails a guard throws
    /// `Failure.nonStandardCell`. Non-finite numbers, occupancy out of range and atoms too close are left to
    /// `CrystalStructure.validationIssues`.
    public static func structure(from doc: MaterialsProjectDocument) throws -> ImportedStructure {
        guard let served = doc.structure else { throw Failure.malformed("\(doc.materialID) has no structure block.") }
        guard let lattice = served.lattice else { throw Failure.malformed("\(doc.materialID)'s structure has no lattice.") }
        guard let rawSites = served.sites, !rawSites.isEmpty else {
            throw Failure.malformed("\(doc.materialID)'s structure lists no atomic sites.")
        }
        var atomSites: [AtomSite] = []
        atomSites.reserveCapacity(rawSites.count)
        for (index, site) in rawSites.enumerated() {
            guard let abc = site.abc, abc.count == 3 else {
                throw Failure.malformed("Site \(index) of \(doc.materialID) has no usable fractional coordinates.")
            }
            guard let species = site.species, !species.isEmpty else {
                throw Failure.malformed("Site \(index) of \(doc.materialID) lists no species.")
            }
            for sp in species {
                let symbol = Element.normalizedSymbol(fromCIF: sp.element)
                guard let z = Element.atomicNumber(forSymbol: symbol) else {
                    throw Failure.malformed(
                        "Site \(index) of \(doc.materialID) has an unrecognized element symbol \"\(sp.element)\".")
                }
                atomSites.append(AtomSite(z: z, fractional: SIMD3(abc[0], abc[1], abc[2]), occupancy: sp.occu,
                                          label: site.label))
            }
        }
        let rows: [SIMD3<Double>]
        if let matrix = lattice.matrix, matrix.count == 3, matrix.allSatisfy({ $0.count == 3 }) {
            rows = matrix.map { SIMD3($0[0], $0[1], $0[2]) }
        } else {
            rows = UnitCell(a: lattice.a, b: lattice.b, c: lattice.c,
                            alphaDeg: lattice.alpha, betaDeg: lattice.beta, gammaDeg: lattice.gamma).realLattice
        }
        let result: CellStandardisation.Result
        do {
            result = try CellStandardisation.standardise(
                latticeRows: rows, sites: atomSites, symbol: doc.symmetry?.symbol, crystalSystem: doc.symmetry?.crystalSystem)
        } catch let failure as CellStandardisation.Failure {
            throw Failure.nonStandardCell("\(doc.materialID): \(failure.message) Import a conventional CIF for this phase instead.")
        }
        let cell: UnitCell
        switch result.outcome {
        case .asServed:
            cell = UnitCell(a: lattice.a, b: lattice.b, c: lattice.c,
                            alphaDeg: lattice.alpha, betaDeg: lattice.beta, gammaDeg: lattice.gamma)
        case .conventional:
            let r = result.latticeRows
            func angle(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double {
                acos(min(1, max(-1, dot(u, v) / (length(u) * length(v))))) * 180 / .pi
            }
            cell = UnitCell(a: length(r[0]), b: length(r[1]), c: length(r[2]),
                            alphaDeg: angle(r[1], r[2]), betaDeg: angle(r[0], r[2]), gammaDeg: angle(r[0], r[1]))
        }
        // MP coordinates are full precision (half-step 0): a trigonal cell in hexagonal axes is not read as 6/mmm, nor Pa-3 as m-3m.
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: result.sites, name: doc.materialID)
        return ImportedStructure(structure: structure, cell: result.outcome)
    }

    /// Grades a fetched document against what the user expected.
    ///
    /// The formula comparison is order-independent by literal element counts ("LiAl2Cu" equals "Al2CuLi", "Al2Cu" does
    /// not equal "AlCu": counts are never reduced to a smallest ratio). A formula match alone is not `.matches` when a
    /// space-group number is also given: both must agree (mp-1185307 was misidentified by formula alone, 2026-09-17).
    public static func check(_ doc: MaterialsProjectDocument, against expectation: PhaseExpectation) -> PhaseExpectationVerdict {
        guard expectation.formula != nil || expectation.spaceGroupNumber != nil else { return .unchecked }
        if let expectedFormula = expectation.formula {
            guard let expectedCounts = parseFormula(expectedFormula) else {
                return .mismatch(reason: "Could not parse expected formula \"\(expectedFormula)\".")
            }
            guard let docFormula = doc.formulaPretty, let docCounts = parseFormula(docFormula) else {
                return .mismatch(reason: "\(doc.materialID) has no usable formula to compare against \"\(expectedFormula)\".")
            }
            guard docCounts == expectedCounts else {
                return .mismatch(reason: "\(doc.materialID) is \(docFormula); expected \(expectedFormula).")
            }
        }
        if let expectedSpaceGroup = expectation.spaceGroupNumber {
            let actual = doc.symmetry?.number
            guard actual == expectedSpaceGroup else {
                let label = doc.formulaPretty ?? doc.materialID
                let description: String
                if let actual {
                    description = doc.symmetry?.symbol.map { "\($0) (\(actual))" } ?? "space group \(actual)"
                } else {
                    description = "an unrecorded space group"
                }
                return .mismatch(reason: "\(doc.materialID) is \(label) in \(description); expected space group \(expectedSpaceGroup).")
            }
        }
        return .matches
    }

    /// Element symbol plus optional integer count, no parentheses ("LiAl2Cu" gives Li 1, Al 2, Cu 1); nil on anything else
    /// (a lowercase-leading token, punctuation, an empty string, a zero or negative count).
    public static func parseFormula(_ text: String) -> [String: Int]? {
        var result: [String: Int] = [:]
        let chars = Array(text)
        guard !chars.isEmpty else { return nil }
        var i = 0
        while i < chars.count {
            guard chars[i].isASCII, chars[i].isUppercase, chars[i].isLetter else { return nil }
            var symbol = String(chars[i])
            i += 1
            while i < chars.count, chars[i].isASCII, chars[i].isLowercase, chars[i].isLetter {
                symbol.append(chars[i])
                i += 1
            }
            var digits = ""
            while i < chars.count, chars[i].isASCII, chars[i].isNumber {
                digits.append(chars[i])
                i += 1
            }
            let count = digits.isEmpty ? 1 : Int(digits)
            guard let count, count > 0 else { return nil }
            result[symbol, default: 0] += count
        }
        return result
    }
}
