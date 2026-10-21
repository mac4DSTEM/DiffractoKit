// swift-tools-version: 6.2
// DiffractoKit benchmark: a SEPARATE package that depends on DiffractoKit by path (like tools/acom-benchmark). Prints JSON lines,
// never asserts. Run with run.sh, under ../../../tools/memcap.sh, one benchmark at a time machine-wide.
import Foundation
import PackageDescription

let kit = URL(fileURLWithPath: Context.packageDirectory).deletingLastPathComponent().deletingLastPathComponent()
    .lastPathComponent

let package = Package(
    name: "DiffractoBench",
    platforms: [.macOS("15.0")],
    dependencies: [.package(path: "../..")],
    targets: [
        .executableTarget(name: "diffracto-bench", dependencies: [
            .product(name: "DiffractoCore", package: kit),
            .product(name: "DiffractoScattering", package: kit),
            .product(name: "DiffractoGeometry", package: kit),
            .product(name: "DiffractoOrientation", package: kit),
        ]),
    ],
    swiftLanguageModes: [.v6]
)
