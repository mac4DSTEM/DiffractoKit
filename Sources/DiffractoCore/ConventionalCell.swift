//
//  ConventionalCell.swift
//  Role: the conventional (possibly centred) cell of a lattice, its Bravais lattice and the integer transformation to
//        it (registration docs/registration-slice-1b.md section 1). The lattice is first Niggli-reduced
//        (`niggliReduced`). On the reduced basis every lattice symmetry is an integer matrix with entries in
//        {-1, 0, 1}, so the lattice's point group is found exactly by testing those 3^9 matrices against the metric
//        (RᵀGR = G within a stated relative tolerance); that group IS the lattice character the International Tables
//        (A, section 9.2) tabulate as 44 reduced-form conditions: its order names the crystal family (2 aP, 4 m, 8 o,
//        16 t, 12 hR, 24 hP, 48 c) and the 44 characters differ only in which centring the family's conventional
//        cell then needs. The conventional axes are the symmetry axes of that group (the shortest lattice vector along
//        each), and the centring is read off the lattice points inside the cell. No table is copied.
//
//        The setting (the choices the International Tables leave open) is fixed here: monoclinic unique axis b with
//        C centring and β ≥ 90°; orthorhombic centring on the (a, b) face; tetragonal and hexagonal unique axis c;
//        hexagonal axes for rhombohedral lattices, obverse; right-handed. Among equivalent choices the one whose
//        transformation is closest to the identity is taken. A cell already in this setting maps to itself
//        (`transformation == .identity`), but a cell that is a conventional cell of its lattice in ANOTHER setting is
//        re-labelled into this one: a monoclinic P cell (index 1) goes to unique axis b with a and c the shortest pair and
//        β ≥ 90°, so a c-unique or a-unique monoclinic cell, or one whose a and c are not that pair, comes back with a
//        permuted or sign-changed transformation (Gate B on 4876b18, F2). Read `transformation` and `setting`, not the
//        family name, to know whether the cell moved.
//
//        Tolerances (Gate B on 4876b18, F3): the lattice symmetry here is decided by the relative `tolerance` alone
//        (|RᵀGR − G|ᵢⱼ ≤ tolerance · √(Gᵢᵢ Gⱼⱼ), default 1e-9). `detectingSymmetry(metricLengthTolerance:)` is passed on as
//        this tolerance; `metricAngleToleranceDeg` is not used by the conventional-cell route to find the cell or its
//        Bravais lattice. A primitive body-centred cubic cell written with the angle 109.4712° (the exact one is
//        109.47122063…°) is therefore a rhombohedral hR lattice at the default and reads Im-3m (cI) only from
//        `metricLengthTolerance` ≈ 1e-5 up (spglib's default symprec of 1e-5 does).
//

import Foundation
import simd

/// The fourteen Bravais lattices.
public enum BravaisLattice: String, Sendable, Equatable, CaseIterable, CustomStringConvertible {
    case aP, mP, mC, oP, oC, oI, oF, tP, tI, hR, hP, cP, cI, cF

    public var description: String { rawValue }

    /// P, C, I, F or R.
    public var centring: Character { rawValue.last! }

    /// The order of the lattice's point group (the holohedry), proper and improper rotations.
    public var holohedryOrder: Int {
        switch self {
        case .aP: return 2
        case .mP, .mC: return 4
        case .oP, .oC, .oI, .oF: return 8
        case .tP, .tI: return 16
        case .hR: return 12
        case .hP: return 24
        case .cP, .cI, .cF: return 48
        }
    }

    /// How many lattice points the conventional cell holds: 1 (P), 2 (C, I), 3 (R in hexagonal axes), 4 (F).
    public var centringMultiplicity: Int {
        switch centring {
        case "P": return 1
        case "C", "I": return 2
        case "R": return 3
        default: return 4
        }
    }
}

/// The conventional cell of a lattice.
public struct ConventionalCell: Sendable, Equatable {
    /// The conventional cell, its parameters idealised to the family's constraints (the six parameters of the input
    /// when `transformation` is the identity).
    public let cell: UnitCell
    public let bravais: BravaisLattice
    /// P: the conventional basis vectors are the COLUMNS of P in terms of the input basis, x_input = P x_conventional
    /// (determinant +multiplicity, so the cell is right-handed). The reflections follow h_conventional = h_input P.
    public let transformation: IntegerMatrix
    public let reduced: ReducedCell
    /// The lattice points inside the conventional cell, in its fractional coordinates, the origin first (one for P).
    public let centringTranslations: [SIMD3<Double>]
    /// The unique axis of a monoclinic (b), tetragonal or hexagonal (c) cell; nil otherwise.
    public let uniqueAxis: UnitCell.Axis?
    /// A one-line statement of the setting the cell is in.
    public let setting: String

    public var centringMultiplicity: Int { transformation.determinant }

    /// Whether the Niggli reduction under this cell reached its conditions (`reduced.converged`). False for an invalid
    /// cell (`UnitCell.isValid`: this includes a cell flatter than V < 1e-6·a·b·c, which is where a shear of about 10^6
    /// lattice units lands), for a shear beyond 2⁴⁰ units, for a reduction that did not stop within its 2000 steps, and
    /// for a cell whose axes differ so much (a billion to one and more) that the integer matrices of the route overflow, and
    /// for a transformation whose determinant is beyond 64 (which a very loose `tolerance` can produce).
    /// Then the conventional route has no reduced basis to work on: `bravais` is `aP` and `transformation` the identity,
    /// which is NOT a finding about the lattice, `setting` says "not reducible", and
    /// `detectingSymmetry(conventional: true)` falls back to the detection in the file's own cell (what
    /// `conventional: false` does) and carries this flag in `ConventionalSetting.reductionConverged`. Read it before
    /// trusting an `aP`. (A shear by any smaller number of lattice units reduces in a few dozen steps; the six
    /// parameters of a cell sheared by k carry a relative rounding of about 1e-16·k², so the default symmetry
    /// tolerance of 1e-9 reads a shear above about 1200 units as the lower family it is within rounding: pass a
    /// larger `tolerance`.)
    public var reductionConverged: Bool { reduced.converged }
}

extension UnitCell {
    /// The default relative tolerance of the lattice-symmetry test (the metric family's 1e-9).
    public static let conventionalCellTolerance = 1e-9

    /// The conventional cell of this lattice (the cell is taken as a primitive lattice) at a relative `tolerance` on the
    /// metric (|RᵀGR − G|ᵢⱼ ≤ tolerance · √(Gᵢᵢ Gⱼⱼ)); this is the only tolerance of the route (the file header, F3: a bcc
    /// primitive cell with 109.4712° reads hR below a tolerance of about 1e-5). An invalid cell, or one whose reduction does
    /// not converge, returns itself as `aP` with `reductionConverged == false`. A cell that is already conventional in
    /// another setting is re-labelled into this one (the file header, F2). A `tolerance` that is NaN, negative or infinite
    /// is read as 0 (the metric must be preserved exactly), never as "every integer matrix is a symmetry"; at 0 the rounding
    /// of an exact cubic cell's parameters can already read it as hR, so pass a real tolerance.
    public func conventionalCell(tolerance: Double = UnitCell.conventionalCellTolerance) -> ConventionalCell {
        LatticeConstruction.conventionalCell(of: self, tolerance: UnitCell.readAsToleranceOrZero(tolerance))
    }
}

// MARK: - Construction

enum LatticeConstruction {

    /// Every integer matrix with entries in {-1, 0, 1} and determinant ±1 (computed once).
    static let smallUnimodular: [IntegerMatrix] = {
        var out: [IntegerMatrix] = []
        let values = [-1, 0, 1]
        var index = 0
        for _ in 0..<19683 {
            var n = index
            var e = [Int](repeating: 0, count: 9)
            for k in 0..<9 { e[k] = values[n % 3]; n /= 3 }
            index += 1
            let m = IntegerMatrix(SIMD3(e[0], e[1], e[2]), SIMD3(e[3], e[4], e[5]), SIMD3(e[6], e[7], e[8]))
            if abs(m.determinant) == 1 { out.append(m) }
        }
        return out
    }()

    /// The integer matrices R with Rᵀ G R = G within `tolerance` (relative): the lattice's point group in the basis of G.
    static func automorphisms(of g: [SIMD3<Double>], tolerance: Double) -> [IntegerMatrix] {
        let scale = [g[0].x.squareRoot(), g[1].y.squareRoot(), g[2].z.squareRoot()]
        var out: [IntegerMatrix] = []
        for r in smallUnimodular {
            let c = (0..<3).map { r.column($0) }
            var ok = true
            search: for i in 0..<3 {
                for j in i..<3 {
                    let h = UnitCell.bilinear(g, c[i], c[j])
                    if abs(h - g[i][j]) > tolerance * scale[i] * scale[j] { ok = false; break search }
                }
            }
            if ok { out.append(r) }
        }
        return out
    }

    static func gcd(_ a: Int, _ b: Int) -> Int { b == 0 ? abs(a) : gcd(b, a % b) }

    static func primitive(_ v: SIMD3<Int>) -> SIMD3<Int> {
        let g = gcd(gcd(v.x, v.y), v.z)
        return g > 1 ? v / g : v
    }

    /// The primitive integer vector along the axis of a proper rotation.
    static func axis(of r: IntegerMatrix) -> SIMD3<Int>? {
        let n = IntegerMatrix(rows: [r.rows[0] &- SIMD3(1, 0, 0), r.rows[1] &- SIMD3(0, 1, 0), r.rows[2] &- SIMD3(0, 0, 1)])
        func cross(_ a: SIMD3<Int>, _ b: SIMD3<Int>) -> SIMD3<Int> {
            SIMD3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
        }
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            let w = cross(n.rows[i], n.rows[j])
            if w != SIMD3(0, 0, 0) {
                var p = primitive(w)
                if let first = [p.x, p.y, p.z].first(where: { $0 != 0 }), first < 0 { p = 0 &- p }
                return p
            }
        }
        return nil
    }

    /// ℓ with {x : ℓ·x = 0} the lattice plane perpendicular to the axis of `r` (order `n`): the rows of Σ Rᵏ are
    /// multiples of ℓ.
    static func planeFunctional(of r: IntegerMatrix, order n: Int) -> SIMD3<Int>? {
        var power = IntegerMatrix.identity
        var sum = [SIMD3<Int>](repeating: SIMD3(0, 0, 0), count: 3)
        for _ in 0..<n {
            for k in 0..<3 { sum[k] &+= power.rows[k] }
            power = power.multiplied(by: r)
        }
        guard let row = sum.first(where: { $0 != SIMD3(0, 0, 0) }) else { return nil }
        var p = primitive(row)
        if let first = [p.x, p.y, p.z].first(where: { $0 != 0 }), first < 0 { p = 0 &- p }
        return p
    }

    /// Two integer vectors that span {x ∈ ℤ³ : ℓ·x = 0} (ℓ primitive).
    static func kernelBasis(of ell: SIMD3<Int>) -> (SIMD3<Int>, SIMD3<Int>) {
        var v = [ell.x, ell.y, ell.z]
        var u: [SIMD3<Int>] = [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)]
        while v.filter({ $0 != 0 }).count > 1 {
            let p = (0..<3).filter { v[$0] != 0 }.min { abs(v[$0]) < abs(v[$1]) }!
            for q in 0..<3 where q != p && v[q] != 0 {
                let k = v[q] / v[p]
                v[q] -= k * v[p]
                u[q] = u[q] &- k &* u[p]
            }
        }
        let keep = (0..<3).filter { v[$0] == 0 }
        return (u[keep[0]], u[keep[1]])
    }

    /// The Gauss-reduced basis of the plane {ℓ·x = 0} under the metric `g`: the shortest vector first.
    static func reducedPlane(_ ell: SIMD3<Int>, _ g: [SIMD3<Double>]) throws -> (SIMD3<Int>, SIMD3<Int>) {
        var (u1, u2) = kernelBasis(of: ell)
        func q(_ x: SIMD3<Int>) -> Double { UnitCell.bilinear(g, x, x) }
        for _ in 0..<100 {
            if q(u2) < q(u1) { swap(&u1, &u2) }
            let mu = (UnitCell.bilinear(g, u1, u2) / q(u1)).rounded()
            if mu == 0 { break }
            // a multiple that is not a number or not an Int (a plane of a billion-to-one cell) is an overflow
            guard abs(mu) < 4.0e18, let next = UnitCell.sheared(u2, by: Int(mu), along: u1) else { throw LatticeOverflow() }
            u2 = next
            if q(u2) >= q(u1) * (1 - 1e-12) { break }
        }
        if q(u2) < q(u1) { swap(&u1, &u2) }
        return (u1, u2)
    }

    /// The distinct lattice points of the conventional cell (numerators over |det P|, each in 0 ..< |det P|).
    static func pointClasses(_ p: IntegerMatrix) throws -> (denominator: Int, classes: Set<SIMD3<Int>>) {
        let d = abs(p.determinant)
        // a conventional cell holds 1 to 4 lattice points (the loop below is d³ long): a determinant beyond 64, or an
        // adjugate that overflows, is not a cell
        guard d <= 64, let adj = p.adjugateChecked else { throw LatticeOverflow() }
        let sign = p.determinant > 0 ? 1 : -1
        var out = Set<SIMD3<Int>>()
        for i in 0..<d { for j in 0..<d { for k in 0..<d {
            let v = adj.times(SIMD3(i, j, k)) &* sign
            out.insert(SIMD3(((v.x % d) + d) % d, ((v.y % d) + d) % d, ((v.z % d) + d) % d))
        } } }
        return (d, out)
    }

    /// Thrown inside the route when an integer product, an adjugate or a plane reduction leaves `Int` (a cell whose axes
    /// differ by a billion to one): the route then answers "not reducible", never traps.
    struct LatticeOverflow: Error {}

    static func fallback(_ cell: UnitCell) -> ConventionalCell {
        let reduced = cell.niggliReduced()
        return ConventionalCell(cell: reduced.cell, bravais: .aP, transformation: reduced.transformation, reduced: reduced,
                                centringTranslations: [SIMD3(0, 0, 0)], uniqueAxis: nil,
                                setting: "aP: the Niggli-reduced cell")
    }

    static func conventionalCell(of input: UnitCell, tolerance: Double) -> ConventionalCell {
        let reduced = input.niggliReduced()
        guard reduced.converged, input.isValid else { return notReducible(input, reduced) }
        do {
            return try build(of: input, reduced: reduced, tolerance: tolerance)
        } catch {
            // an integer overflow: the cell is flagged as not reduced (reductionConverged == false)
            let flagged = ReducedCell(cell: reduced.cell, transformation: reduced.transformation, epsilon: reduced.epsilon,
                                      converged: false, steps: reduced.steps)
            return notReducible(input, flagged)
        }
    }

    static func notReducible(_ input: UnitCell, _ reduced: ReducedCell) -> ConventionalCell {
        ConventionalCell(cell: input, bravais: .aP, transformation: .identity, reduced: reduced,
                         centringTranslations: [SIMD3(0, 0, 0)], uniqueAxis: nil,
                         setting: "aP: the cell as given (not reducible)")
    }

    static func build(of input: UnitCell, reduced: ReducedCell, tolerance: Double) throws -> ConventionalCell {
        let pRed = reduced.transformation
        let g = UnitCell.gram(of: pRed, input.metricTensor)         // the reduced metric, exact from the input's
        let group = automorphisms(of: g, tolerance: tolerance)
        let proper = group.filter { $0.determinant == 1 }
        func q(_ x: SIMD3<Int>) -> Double { UnitCell.bilinear(g, x, x) }

        var cols: [SIMD3<Int>]
        var unique: UnitCell.Axis?
        var description: String
        switch group.count {
        case 48:
            var axes: [SIMD3<Int>] = []
            for r in proper where r.order == 4 {
                if let w = axis(of: r), !axes.contains(w) { axes.append(w) }
            }
            guard axes.count == 3 else { return fallback(input) }
            cols = axes
            description = "cubic axes along the three 4-fold axes"
        case 24:
            guard let six = proper.first(where: { $0.order == 6 }), let w = axis(of: six),
                  let ell = planeFunctional(of: six, order: 6) else { return fallback(input) }
            let (p, _) = try reducedPlane(ell, g)
            let three = six.multiplied(by: six)
            cols = [p, three.times(p), w]
            unique = .c
            description = "hexagonal axes, unique axis c, γ = 120°"
        case 12:
            guard let three = proper.first(where: { $0.order == 3 }), let w = axis(of: three),
                  let ell = planeFunctional(of: three, order: 3) else { return fallback(input) }
            let (p, _) = try reducedPlane(ell, g)
            cols = [p, three.times(p), w]
            unique = .c
            description = "hexagonal axes of a rhombohedral lattice, obverse, γ = 120°"
        case 16:
            guard let four = proper.first(where: { $0.order == 4 }), let w = axis(of: four),
                  let ell = planeFunctional(of: four, order: 4) else { return fallback(input) }
            let (p, _) = try reducedPlane(ell, g)
            cols = [p, four.times(p), w]
            unique = .c
            description = "tetragonal, unique axis c"
        case 8:
            let axes = proper.filter { $0.order == 2 }.compactMap { axis(of: $0) }
            guard axes.count == 3 else { return fallback(input) }
            cols = axes.sorted { q($0) < q($1) }
            description = "orthorhombic along the three 2-fold axes"
        case 4:
            guard let two = proper.first(where: { $0.order == 2 }), let w = axis(of: two),
                  let ell = planeFunctional(of: two, order: 2) else { return fallback(input) }
            let (e1, e2) = try reducedPlane(ell, g)
            unique = .b
            let index = abs(IntegerMatrix(columns: [e1, e2, w]).determinant)
            if index == 1 {
                var c = e2
                if UnitCell.bilinear(g, e1, c) > 0 { c = 0 &- c }
                cols = [e1, w, c]
                description = "monoclinic P, unique axis b, a and c the shortest pair with β ≥ 90°"
            } else if index == 2 {
                // the centring class: λ = [e1 e2 w]⁻¹ n mod 1 for a lattice vector not in the plane+axis sublattice
                let basis = IntegerMatrix(columns: [e1, e2, w])
                guard let adj = basis.adjugateChecked else { throw LatticeOverflow() }
                let sign = basis.determinant > 0 ? 1 : -1
                var parity = SIMD2<Int>(0, 0)
                for n in [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)] {
                    let v = adj.times(n) &* sign
                    let m = SIMD3(((v.x % 2) + 2) % 2, ((v.y % 2) + 2) % 2, ((v.z % 2) + 2) % 2)
                    if m.z == 1 { parity = SIMD2(m.x, m.y); break }
                }
                // a ≡ 2u (mod 2Λ), primitive, shortest
                var best: (SIMD3<Int>, Double)?
                var bestCoefficients = SIMD2<Int>(0, 0)
                for m1 in -3...3 { for m2 in -3...3 where abs(gcd(m1, m2)) == 1
                    && ((m1 % 2) + 2) % 2 == parity.x && ((m2 % 2) + 2) % 2 == parity.y {
                    let a = m1 &* e1 &+ m2 &* e2
                    let len = q(a)
                    if best == nil || len < best!.1 * (1 - 1e-9) { best = (a, len); bestCoefficients = SIMD2(m1, m2) }
                } }
                guard let (a, _) = best else { return fallback(input) }
                var c: SIMD3<Int>?
                var cLen = Double.infinity
                for k1 in -3...3 { for k2 in -3...3
                    where abs(bestCoefficients.x * k2 - bestCoefficients.y * k1) == 1 {
                    let cand = k1 &* e1 &+ k2 &* e2
                    if q(cand) < cLen * (1 - 1e-9) { c = cand; cLen = q(cand) }
                } }
                guard var cc = c else { return fallback(input) }
                if UnitCell.bilinear(g, a, cc) > 0 { cc = 0 &- cc }
                cols = [a, w, cc]
                description = "monoclinic C, unique axis b, C centring, a the shortest lattice vector with the centring along it"
            } else { return fallback(input) }
        default:
            return fallback(input)
        }

        var pc = IntegerMatrix(columns: cols)
        // the family's conventional cell can still be P, centred, or (orthorhombic) differently centred
        var (d, classes) = try pointClasses(pc)
        var bravais: BravaisLattice
        switch group.count {
        case 48: bravais = d == 1 ? .cP : (d == 2 ? .cI : .cF)
        case 24: bravais = .hP
        case 12: bravais = .hR
        case 16: bravais = d == 1 ? .tP : .tI
        case 8:
            if d == 1 { bravais = .oP } else if d == 4 { bravais = .oF } else if classes.contains(SIMD3(1, 1, 1)) { bravais = .oI } else {
                bravais = .oC
                // bring the centred face to (a, b): cyclic permutation of the axes
                let nonzero = classes.first { $0 != SIMD3(0, 0, 0) }!
                let z = [nonzero.x, nonzero.y, nonzero.z].firstIndex(of: 0)!
                let order = z == 0 ? [1, 2, 0] : (z == 1 ? [2, 0, 1] : [0, 1, 2])
                cols = order.map { cols[$0] }
                if q(cols[0]) > q(cols[1]) { cols.swapAt(0, 1) }
                pc = IntegerMatrix(columns: cols)
            }
        case 4: bravais = d == 1 ? .mP : .mC
        default: bravais = .aP
        }
        if group.count == 2 { cols = [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)]; pc = IntegerMatrix(columns: cols) }
        if group.count == 12 {
            guard pc.determinant != 0, abs(pc.determinant) == 3 else { return fallback(input) }
            // reverse setting → obverse: (a, b, c) → (−a, −b, c)
            if !(try pointClasses(pc)).classes.contains(SIMD3(2, 1, 1)) {
                cols[0] = 0 &- cols[0]; cols[1] = 0 &- cols[1]; pc = IntegerMatrix(columns: cols)
            }
        }
        if pc.determinant < 0 { pc = pc.negated }
        (d, classes) = try pointClasses(pc)

        // Among the equivalent conventional bases (the cell's own symmetries that keep the lattice, and for orthorhombic
        // every signed permutation that keeps the centring), the one closest to the identity.
        let gc = UnitCell.gram(of: pc, g)
        var candidates = automorphisms(of: gc, tolerance: tolerance).filter { $0.determinant == 1 }
        if bravais.rawValue.hasPrefix("o") {
            let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
            for perm in perms { for s0 in [-1, 1] { for s1 in [-1, 1] { for s2 in [-1, 1] {
                let signs = [s0, s1, s2]
                var rows: [SIMD3<Int>] = []
                for i in 0..<3 { var row = SIMD3(0, 0, 0); row[perm[i]] = signs[i]; rows.append(row) }
                let m = IntegerMatrix(rows: rows)
                if m.determinant == 1 && !candidates.contains(m) { candidates.append(m) }
            } } } }
        }
        var bestMatrix = pc
        var bestScore = Int.max
        var bestKey: [Int] = []
        for candidate in candidates {
            guard let trial = pc.multipliedChecked(by: candidate) else { throw LatticeOverflow() }
            let trialClasses = try pointClasses(trial)
            guard trialClasses.denominator == d, trialClasses.classes == classes else { continue }
            guard let total = pRed.multipliedChecked(by: trial) else { throw LatticeOverflow() }
            var score = 0
            var key: [Int] = []
            for i in 0..<3 { for j in 0..<3 { score += abs(total[i, j] - (i == j ? 1 : 0)); key.append(total[i, j]) } }
            if score < bestScore || (score == bestScore && key.lexicographicallyPrecedes(bestKey)) {
                bestScore = score; bestKey = key; bestMatrix = trial
            }
        }
        pc = bestMatrix
        guard let total = pRed.multipliedChecked(by: pc), total.adjugateChecked != nil else { throw LatticeOverflow() }
        let gTotal = UnitCell.gram(of: total, input.metricTensor)
        let cell: UnitCell
        if total == .identity { cell = input } else { cell = idealised(gTotal, bravais) }
        let translations = classes.sorted { l, r in
            let lz = l == SIMD3(0, 0, 0), rz = r == SIMD3(0, 0, 0)
            if lz != rz { return lz }
            return [l.x, l.y, l.z].lexicographicallyPrecedes([r.x, r.y, r.z])
        }.map { SIMD3(Double($0.x), Double($0.y), Double($0.z)) / Double(d) }
        let fullDescription = "\(bravais.rawValue): \(description)"
        return ConventionalCell(cell: cell, bravais: bravais, transformation: total, reduced: reduced,
                                centringTranslations: translations, uniqueAxis: unique, setting: fullDescription)
    }

    /// The cell of metric `g` with the family's constraints imposed exactly.
    static func idealised(_ g: [SIMD3<Double>], _ bravais: BravaisLattice) -> UnitCell {
        let l = [g[0].x.squareRoot(), g[1].y.squareRoot(), g[2].z.squareRoot()]
        switch bravais {
        case .cP, .cI, .cF: return .cubic(a: l[0])
        case .tP, .tI: return .tetragonal(a: l[0], c: l[2])
        case .hP, .hR: return .hexagonal(a: l[0], c: l[2])
        case .oP, .oC, .oI, .oF: return .orthorhombic(a: l[0], b: l[1], c: l[2])
        case .mP, .mC:
            let beta = acos(min(1, max(-1, g[0].z / (l[0] * l[2])))) * 180 / .pi
            return .monoclinic(a: l[0], b: l[1], c: l[2], betaDeg: beta)
        case .aP: return UnitCell.cell(fromMetric: g)
        }
    }
}
