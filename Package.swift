// swift-tools-version: 6.2

import PackageDescription

let package = Package(
    name: "MayorFEM",
    platforms: [
        .macOS(.v13),
    ],
    products: [
        .library(name: "MayorFEM", targets: ["MayorFEM"]),
        .executable(name: "mayor-fem", targets: ["MayorFEMCLI"]),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-testing.git", from: "0.12.0"),
    ],
    targets: [
        .target(
            name: "MayorFEM",
            resources: [
                .process("Resources"),
            ]
        ),
        .executableTarget(
            name: "MayorFEMCLI",
            dependencies: ["MayorFEM"]
        ),
        .testTarget(
            name: "MayorFEMTests",
            dependencies: [
                "MayorFEM",
                .product(name: "Testing", package: "swift-testing"),
            ]
        ),
    ]
)
