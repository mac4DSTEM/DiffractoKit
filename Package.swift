// swift-tools-version: 6.2
// DiffractoKit — the macEM suite's crystallography and diffraction library: cells, symmetry, CIF import, d-spacings and
// reflection orbits (DiffractoCore, slice 1); scattering and structure factors (DiffractoScattering, slice 2); pattern
// geometry, calibration, camera length, zone axes and indexing (DiffractoGeometry); orientation results and Euler
// conventions (DiffractoOrientation). The suite ROADMAP.md "Kits" row is the direction; docs/registration-slice-1.md the
// first registration. Ported from mac4DSTEM 7448401b and macTEM a58c84c (both GPL-3.0), checked against pymatgen (MIT) as
// an oracle; no pymatgen code is copied.
//
// Swift 6 language mode, nonisolated by default (no `.defaultIsolation`), the three upcoming features the suite's
// packages enable, no unsafe flags (a versioned dependency cannot carry them): a consumer builds the kit optimised
// itself, and tools/test.sh runs this package's tests with `-c release`. One product per module. What one module needs
// from another is `package`, never `public`; what an app reads is `public` and is guarded by tools/consumer-check.
import PackageDescription

let settings: [SwiftSetting] = [
    .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
    .enableUpcomingFeature("InferIsolatedConformances"),
    .enableUpcomingFeature("MemberImportVisibility"),
]

let package = Package(
    name: "DiffractoKit",
    platforms: [.macOS("27.0")],
    products: [
        .library(name: "DiffractoCore", targets: ["DiffractoCore"]),
        .library(name: "DiffractoScattering", targets: ["DiffractoScattering"]),
        .library(name: "DiffractoGeometry", targets: ["DiffractoGeometry"]),
        .library(name: "DiffractoOrientation", targets: ["DiffractoOrientation"]),
    ],
    targets: [
        // Slice 1 (docs/registration-slice-1.md): cells of the seven systems, symmetry operations and Laue groups, CIF
        // import, the reflection list with systematic absences, orbits and families. Foundation and simd only.
        .target(name: "DiffractoCore", swiftSettings: settings),
        // Slice 2 (next registration): electron scattering factors (Lobato–Van Dyck, py4DSTEM's table), structure factors
        // and kinematic intensities, against pymatgen and py4DSTEM's tables. Empty until registered.
        .target(name: "DiffractoScattering", dependencies: ["DiffractoCore"], swiftSettings: settings),
        // Later: patterns, ring and ellipse calibration, camera length, zone axes, indexing. Empty until registered.
        .target(name: "DiffractoGeometry", dependencies: ["DiffractoCore"], swiftSettings: settings),
        // Later: orientation results, symmetry reduction, Euler conventions against pyxem / orix. Empty until registered.
        .target(name: "DiffractoOrientation", dependencies: ["DiffractoCore"], swiftSettings: settings),
        .testTarget(name: "DiffractoCoreTests", dependencies: ["DiffractoCore"], exclude: ["Fixtures"],
                    swiftSettings: settings),
        .testTarget(name: "DiffractoScatteringTests", dependencies: ["DiffractoScattering"], swiftSettings: settings),
        .testTarget(name: "DiffractoGeometryTests", dependencies: ["DiffractoGeometry"], swiftSettings: settings),
        .testTarget(name: "DiffractoOrientationTests", dependencies: ["DiffractoOrientation"], swiftSettings: settings),
    ],
    swiftLanguageModes: [.v6]
)
