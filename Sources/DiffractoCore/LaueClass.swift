//
//  LaueClass.swift
//  Role: the seven crystal systems and the eleven Laue classes (the point groups a kinematic pattern can tell apart,
//        since Friedel's law adds the inversion), with the 32 → 11 table and the properties the orbit naming needs.
//

import Foundation

public enum CrystalSystem: String, Sendable, Equatable, CaseIterable, CustomStringConvertible {
    case triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic

    public var description: String { rawValue }
}

/// The eleven Laue classes, named by their Hermann–Mauguin symbol.
public enum LaueClass: String, Sendable, Equatable, CaseIterable, CustomStringConvertible {
    case oneBar = "-1"
    case twoOverM = "2/m"
    case mmm = "mmm"
    case fourOverM = "4/m"
    case fourOverMmm = "4/mmm"
    case threeBar = "-3"
    case threeBarM = "-3m"
    case sixOverM = "6/m"
    case sixOverMmm = "6/mmm"
    case mThreeBar = "m-3"
    case mThreeBarM = "m-3m"

    /// The Hermann–Mauguin symbol.
    public var symbol: String { rawValue }

    /// The order of the Laue group (its proper and improper rotations).
    public var order: Int {
        switch self {
        case .oneBar: return 2
        case .twoOverM: return 4
        case .mmm, .fourOverM: return 8
        case .fourOverMmm: return 16
        case .threeBar: return 6
        case .threeBarM, .sixOverM: return 12
        case .sixOverMmm, .mThreeBar: return 24
        case .mThreeBarM: return 48
        }
    }

    public var crystalSystem: CrystalSystem {
        switch self {
        case .oneBar: return .triclinic
        case .twoOverM: return .monoclinic
        case .mmm: return .orthorhombic
        case .fourOverM, .fourOverMmm: return .tetragonal
        case .threeBar, .threeBarM: return .trigonal
        case .sixOverM, .sixOverMmm: return .hexagonal
        case .mThreeBar, .mThreeBarM: return .cubic
        }
    }

    /// Whether orbit names print four Miller–Bravais indices (the hexagonal-axes classes).
    public var usesMillerBravais: Bool {
        switch self {
        case .threeBar, .threeBarM, .sixOverM, .sixOverMmm: return true
        default: return false
        }
    }

    /// The Laue class of a crystallographic point group, by its Hermann–Mauguin symbol as pymatgen and the IUCr write
    /// it (`"32"`, `"-42m"`, `"m-3"`, `"6/mmm"`); spaces and the overbar spellings `-3`, `3̄` are accepted. nil for an
    /// unknown symbol.
    public init?(pointGroupSymbol raw: String) {
        let symbol = raw.replacingOccurrences(of: " ", with: "").replacingOccurrences(of: "\u{0304}", with: "")
            .replacingOccurrences(of: "̄", with: "")
        let table: [String: LaueClass] = [
            "1": .oneBar, "-1": .oneBar,
            "2": .twoOverM, "m": .twoOverM, "2/m": .twoOverM,
            "222": .mmm, "mm2": .mmm, "mmm": .mmm,
            "4": .fourOverM, "-4": .fourOverM, "4/m": .fourOverM,
            "422": .fourOverMmm, "4mm": .fourOverMmm, "-42m": .fourOverMmm, "-4m2": .fourOverMmm,
            "4/mmm": .fourOverMmm,
            "3": .threeBar, "-3": .threeBar,
            "32": .threeBarM, "3m": .threeBarM, "-3m": .threeBarM, "321": .threeBarM, "312": .threeBarM,
            "3m1": .threeBarM, "31m": .threeBarM, "-3m1": .threeBarM, "-31m": .threeBarM,
            "6": .sixOverM, "-6": .sixOverM, "6/m": .sixOverM,
            "622": .sixOverMmm, "6mm": .sixOverMmm, "-6m2": .sixOverMmm, "-62m": .sixOverMmm, "6/mmm": .sixOverMmm,
            "23": .mThreeBar, "m-3": .mThreeBar, "m3": .mThreeBar,
            "432": .mThreeBarM, "-43m": .mThreeBarM, "m-3m": .mThreeBarM, "m3m": .mThreeBarM,
        ]
        guard let found = table[symbol] else { return nil }
        self = found
    }

    public var description: String { rawValue }
}
