// DiffractoKit slice 2, oracle 3 (registration docs/registration-slice-2.md section 2 and 4.4): a scratch tool that links the
// apps' own Crystal.swift and ScatteringFactors.swift (mac4DSTEM 7d6fdbd1, fetched by run.sh with `git show`, compiled
// beside this file) and dumps what `Crystal.reflections(kMax:tolerance: 1e-4)` returns for the apps' presets: the cell, the
// site list in the apps' order, and for every reflection h, k, l and the Double bit patterns of |g|, Re F, Im F and |F|^2.
//
//   apps-parity <out.json>
//
// Written from the registration, not from the kit. No kit code is read or linked. Decimal values use Swift's shortest
// round-trip description, so they parse back to the same Double; the bit patterns are the evidence.

import Foundation

func hex(_ x: Double) -> String {
    let s = String(x.bitPattern, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}

func dec(_ x: Double) -> String {
    precondition(x.isFinite)
    return "\(x)"
}

struct Case {
    let name: String
    let crystal: Crystal
    let kMax: Double
}

let presets: [(String, Crystal)] = [
    ("aluminum", .aluminum), ("gold", .gold), ("nickel", .nickel), ("copper", .copper), ("iron", .iron),
    ("silicon", .silicon), ("magnesium", .magnesium), ("tungstenDisulfide", .tungstenDisulfide),
    ("betaDoublePrime", .betaDoublePrime),
]
var cases: [Case] = presets.map { Case(name: $0.0, crystal: $0.1, kMax: 1.2) }
cases.append(Case(name: "betaDoublePrime", crystal: .betaDoublePrime, kMax: 2.0))

var out = "{\"tolerance\":1e-4,\"unsupportedElements_all_empty\":"
out += cases.allSatisfy { $0.crystal.unsupportedElements.isEmpty } ? "true" : "false"
out += ",\"cases\":[\n"
for (ci, c) in cases.enumerated() {
    let cr = c.crystal
    let refl = cr.reflections(kMax: c.kMax, tolerance: 1e-4)
    out += "{\"preset\":\"\(c.name)\",\"kMax\":\(dec(c.kMax)),\"count\":\(refl.count),"
    out += "\"cell\":{\"a\":\(dec(cr.a)),\"b\":\(dec(cr.b)),\"c\":\(dec(cr.c)),\"alpha\":\(dec(cr.alphaDeg)),\"beta\":\(dec(cr.betaDeg)),\"gamma\":\(dec(cr.gammaDeg))},"
    out += "\"volume\":\(dec(cr.volume)),\"volumeBits\":\"\(hex(cr.volume))\","
    out += "\"sites\":[" + cr.sites.map { "[\($0.z),\(dec($0.fractional.x)),\(dec($0.fractional.y)),\(dec($0.fractional.z)),\(dec($0.occupancy))]" }.joined(separator: ",") + "],"
    out += "\"siteBits\":[" + cr.sites.map { "[\"\(hex($0.fractional.x))\",\"\(hex($0.fractional.y))\",\"\(hex($0.fractional.z))\",\"\(hex($0.occupancy))\"]" }.joined(separator: ",") + "],"
    out += "\"reflections\":[\n"
    var rows: [String] = []
    rows.reserveCapacity(refl.count)
    for r in refl {
        rows.append("{\"h\":\(r.h),\"k\":\(r.k),\"l\":\(r.l),\"g\":\"\(hex(r.gLength))\",\"re\":\"\(hex(r.structureFactorRe))\",\"im\":\"\(hex(r.structureFactorIm))\",\"i\":\"\(hex(r.intensity))\",\"gv\":\(dec(r.gLength)),\"rev\":\(dec(r.structureFactorRe)),\"imv\":\(dec(r.structureFactorIm)),\"iv\":\(dec(r.intensity))}")
    }
    out += rows.joined(separator: ",\n")
    out += "]}" + (ci + 1 < cases.count ? ",\n" : "\n")
}
out += "]}\n"

let path = CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "slice-2-apps.json"
try out.write(toFile: path, atomically: true, encoding: .utf8)
print("apps-parity: wrote \(path), " + cases.map { "\($0.name)@\($0.kMax)=\($0.crystal.reflections(kMax: $0.kMax, tolerance: 1e-4).count)" }.joined(separator: " "))
