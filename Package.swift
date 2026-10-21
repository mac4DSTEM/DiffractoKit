// swift-tools-version: 6.2
// DiffractoKit — a crystallography and diffraction library: cells, symmetry, CIF import, d-spacings and
// reflection orbits (DiffractoCore); scattering and structure factors (DiffractoScattering); pattern
// geometry, calibration, camera length, zone axes and indexing (DiffractoGeometry); orientation results and Euler
// conventions (DiffractoOrientation). Ported from earlier GPL-3.0 code (commits named in NOTICE), checked against pymatgen (MIT) as a
// reference; no pymatgen code is copied.
//
// Swift 6 language mode, nonisolated by default (no `.defaultIsolation`), the three upcoming features consuming packages
// typically enable, no unsafe flags (a versioned dependency cannot carry them): a consumer builds the kit optimised
// itself, and tools/test.sh runs this package's tests with `-c release`. One product per module. What one module needs
// from another is `package`, never `public`; what a consumer reads is `public` and is guarded by tools/consumer-check.
import PackageDescription

let settings: [SwiftSetting] = [
    .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
    .enableUpcomingFeature("InferIsolatedConformances"),
    .enableUpcomingFeature("MemberImportVisibility"),
]

let package = Package(
    name: "DiffractoKit",
    platforms: [.macOS("15.0"), .iOS("18.0")],
    products: [
        .library(name: "DiffractoCore", targets: ["DiffractoCore"]),
        .library(name: "DiffractoScattering", targets: ["DiffractoScattering"]),
        .library(name: "DiffractoGeometry", targets: ["DiffractoGeometry"]),
        .library(name: "DiffractoOrientation", targets: ["DiffractoOrientation"]),
    ],
    dependencies: [
        // FFT1D for the ACOM matcher's azimuthal ring transforms (DiffractoOrientation only); a range below 1.0.0 (tested against 0.2.0 and 0.4.0); consumers who pin exact tags resolve one ImagoKit.
        .package(url: "https://github.com/mac4DSTEM/ImagoKit", "0.2.0"..<"1.0.0"),
    ],
    targets: [
        // Cells of the seven systems, symmetry operations and Laue groups, CIF
        // import, the reflection list with systematic absences, orbits and families. Foundation, simd and CryptoKit (the CIF bytes hash).
        .target(name: "DiffractoCore", swiftSettings: settings),
        // Electron scattering factors (Lobato–Van Dyck, py4DSTEM's table),
        // structure factors and the kinematic ring list after the |F| cut.
        .target(name: "DiffractoScattering", dependencies: ["DiffractoCore"], swiftSettings: settings),
        // Wavelength, ring calibration and camera length, the indexing decision,
        // zone axes; the ring and ellipse fit.
        .target(name: "DiffractoGeometry", dependencies: ["DiffractoCore", "DiffractoScattering"], swiftSettings: settings),
        // The reference library and the vector matcher of phase mapping, the value types of an orientation result (Euler
        // angles, the zone-axis frame) and the ACOM template matcher. Later: py4DSTEM's SO(3) correlogram and pyxem / orix
        // template matching.
        .target(name: "DiffractoOrientation",
                dependencies: ["DiffractoCore", "DiffractoScattering",
                               .product(name: "ImagoFFT", package: "ImagoKit")],
                swiftSettings: settings),
        // The README's first ten minutes as a program (no product: `swift run diffracto-example`, or `swift build --product
        // diffracto-example`).
        .executableTarget(name: "diffracto-example",
                          dependencies: ["DiffractoCore", "DiffractoScattering", "DiffractoGeometry"],
                          path: "Examples/diffracto-example", swiftSettings: settings),
        .testTarget(name: "DiffractoCoreTests", dependencies: ["DiffractoCore"], exclude: ["Fixtures"],
                    swiftSettings: settings),
        .testTarget(name: "DiffractoScatteringTests", dependencies: ["DiffractoScattering", "DiffractoCore"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "DiffractoGeometryTests", dependencies: ["DiffractoGeometry", "DiffractoScattering", "DiffractoCore"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "DiffractoOrientationTests",
                    dependencies: ["DiffractoOrientation", "DiffractoScattering", "DiffractoCore"],
                    exclude: ["Fixtures"], swiftSettings: settings),
    ],
    swiftLanguageModes: [.v6]
)
