// swift-tools-version: 6.2
// DiffractoKit consumer check: a SEPARATE package that depends on DiffractoKit by path and uses, from outside the
// package, the public API an app reads. `package` access does not cross a package boundary, so this compiles only while
// each symbol it names is `public`. See README.md.
import Foundation
import PackageDescription

// A path dependency's identity is its directory's name (lower-cased): "diffractokit" in the repository, but not in a
// scratch copy; `.product` needs that name, so it is taken from the directory two levels above this manifest.
let kit = URL(fileURLWithPath: Context.packageDirectory).deletingLastPathComponent().deletingLastPathComponent()
    .lastPathComponent

let package = Package(
    name: "ConsumerCheck",
    platforms: [.macOS("15.0")],
    dependencies: [
        .package(path: "../.."),
    ],
    targets: [
        .executableTarget(name: "consumer-check", dependencies: [
            .product(name: "DiffractoCore", package: kit),
            .product(name: "DiffractoScattering", package: kit),
            .product(name: "DiffractoGeometry", package: kit),
            .product(name: "DiffractoOrientation", package: kit),
        ]),
    ],
    swiftLanguageModes: [.v6]
)
