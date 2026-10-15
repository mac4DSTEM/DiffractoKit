//
//  Element.swift
//  Role: chemical symbols 1 ... 118 and the CIF symbol rule (reference commit 7448401b `CIFImport.elementSymbol(from:)`):
//        trailing oxidation-state and isotope markers stripped, case normalised. Which elements have a scattering factor
//        is slice 2's table to say, not this file's.
//

import Foundation

public enum Element {
    /// Symbols by atomic number, 1 ... 118.
    public static let symbols: [String] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
        "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
        "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
        "Lv", "Ts", "Og",
    ]

    public static let maximumAtomicNumber = 118

    private static let atomicNumbers: [String: Int] = {
        var map: [String: Int] = [:]
        for (index, symbol) in symbols.enumerated() { map[symbol] = index + 1 }
        return map
    }()

    /// "Fe" for 26; nil outside 1 ... 118.
    public static func symbol(forAtomicNumber z: Int) -> String? {
        guard z >= 1, z <= symbols.count else { return nil }
        return symbols[z - 1]
    }

    /// 26 for "Fe" (exact symbol); nil for an unknown one. Use `normalizedSymbol(fromCIF:)` first on CIF text.
    public static func atomicNumber(forSymbol symbol: String) -> Int? {
        atomicNumbers[symbol]
    }

    /// Strips trailing oxidation-state and isotope markers and normalises case: `Si1` → `Si`, `Fe3+` → `Fe`,
    /// `O2-` → `O`, `SI` → `Si`.
    public static func normalizedSymbol(fromCIF raw: String) -> String {
        var chars = Array(raw)
        while let last = chars.last, last.isNumber || last == "+" || last == "-" {
            chars.removeLast()
        }
        let stripped = String(chars)
        guard let first = stripped.first else { return stripped }
        let rest = stripped.dropFirst().lowercased()
        return String(first).uppercased() + rest
    }
}
