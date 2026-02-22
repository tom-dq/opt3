import Foundation
import simd

public enum TopologyLiteratureProblems2D {
    public enum LBracketLoadLocation {
        case topRight
        case midRight
    }

    public static func amirLBracket(
        loadLocation: LBracketLoadLocation,
        nx: Int = 24,
        ny: Int = 24,
        order: ElementOrder2D = .linear,
        integrationScheme: IntegrationScheme2D = .full,
        subdivisionLevels: Int = 0,
        endDisplacement: Float = 0.05,
        loadSteps: Int = 10
    ) -> FEMProblem2D {
        var mesh = Mesh2D.rectangularPlate(length: 1.0, height: 1.0, nx: nx, ny: ny, order: order)
        mesh = filteredQuadMesh(mesh: mesh) { centroid in
            !(centroid.x > 0.5 && centroid.y > 0.5)
        }
        if subdivisionLevels > 0 {
            mesh = mesh.subdivided(levels: subdivisionLevels)
        }

        let bounds = meshBounds(mesh)
        let tol = coordinateTolerance(bounds: bounds)
        let xMax = bounds.max.x
        let yMax = bounds.max.y

        let loadedNodes: [Int]
        switch loadLocation {
        case .topRight:
            loadedNodes = boundaryNodes(
                in: mesh,
                where: { point in
                    abs(point.x - xMax) <= tol && point.y >= yMax * 0.82
                },
                fallbackTarget: SIMD2<Float>(xMax, yMax * 0.9),
                fallbackCount: 8
            )
        case .midRight:
            loadedNodes = boundaryNodes(
                in: mesh,
                where: { point in
                    abs(point.x - xMax) <= tol && point.y >= yMax * 0.40 && point.y <= yMax * 0.60
                },
                fallbackTarget: SIMD2<Float>(xMax, yMax * 0.5),
                fallbackCount: 8
            )
        }

        let prescribed = makePrescribedDisplacements(
            mesh: mesh,
            fixedNodePredicate: { point in abs(point.x - bounds.min.x) <= tol },
            loadedNodes: loadedNodes,
            loadedDisplacement: SIMD2<Float>(0, -abs(endDisplacement))
        )

        return FEMProblem2D(
            mesh: mesh,
            material: amirMaterial(),
            thickness: 1.0,
            integrationScheme: integrationScheme,
            prescribedDisplacements: prescribed,
            controls: benchmarkSolverControls(loadSteps: loadSteps)
        )
    }

    public static func amirUBracket(
        nx: Int = 28,
        ny: Int = 20,
        order: ElementOrder2D = .linear,
        integrationScheme: IntegrationScheme2D = .full,
        subdivisionLevels: Int = 0,
        endDisplacement: Float = 0.045,
        loadSteps: Int = 10
    ) -> FEMProblem2D {
        var mesh = Mesh2D.rectangularPlate(length: 1.2, height: 1.0, nx: nx, ny: ny, order: order)
        mesh = filteredQuadMesh(mesh: mesh) { centroid in
            let insideCutoutX = centroid.x >= 0.42 && centroid.x <= 0.78
            let insideCutoutY = centroid.y >= 0.38
            return !(insideCutoutX && insideCutoutY)
        }
        if subdivisionLevels > 0 {
            mesh = mesh.subdivided(levels: subdivisionLevels)
        }

        let bounds = meshBounds(mesh)
        let tol = coordinateTolerance(bounds: bounds)
        let xMax = bounds.max.x
        let yMax = bounds.max.y

        let loadedNodes = boundaryNodes(
            in: mesh,
            where: { point in
                abs(point.x - xMax) <= tol && point.y >= yMax * 0.78
            },
            fallbackTarget: SIMD2<Float>(xMax, yMax * 0.9),
            fallbackCount: 10
        )

        let prescribed = makePrescribedDisplacements(
            mesh: mesh,
            fixedNodePredicate: { point in abs(point.x - bounds.min.x) <= tol },
            loadedNodes: loadedNodes,
            loadedDisplacement: SIMD2<Float>(abs(endDisplacement), 0)
        )

        return FEMProblem2D(
            mesh: mesh,
            material: amirMaterial(),
            thickness: 1.0,
            integrationScheme: integrationScheme,
            prescribedDisplacements: prescribed,
            controls: benchmarkSolverControls(loadSteps: loadSteps)
        )
    }

    public static func liangCantilever(
        loadMagnitude: Float,
        nx: Int = 36,
        ny: Int = 12,
        order: ElementOrder2D = .linear,
        integrationScheme: IntegrationScheme2D = .full,
        subdivisionLevels: Int = 0,
        loadSteps: Int = 10
    ) -> FEMProblem2D {
        let positiveLoad = max(0.5, loadMagnitude)
        let endDisplacement = 2.0e-5 * positiveLoad

        var mesh = Mesh2D.rectangularPlate(length: 0.12, height: 0.03, nx: nx, ny: ny, order: order)
        if subdivisionLevels > 0 {
            mesh = mesh.subdivided(levels: subdivisionLevels)
        }

        let bounds = meshBounds(mesh)
        let tol = coordinateTolerance(bounds: bounds)
        let xMax = bounds.max.x
        let midY = 0.5 * (bounds.min.y + bounds.max.y)
        let halfBand = 0.12 * (bounds.max.y - bounds.min.y)

        let loadedNodes = boundaryNodes(
            in: mesh,
            where: { point in
                abs(point.x - xMax) <= tol && abs(point.y - midY) <= halfBand
            },
            fallbackTarget: SIMD2<Float>(xMax, midY),
            fallbackCount: 8
        )

        let prescribed = makePrescribedDisplacements(
            mesh: mesh,
            fixedNodePredicate: { point in abs(point.x - bounds.min.x) <= tol },
            loadedNodes: loadedNodes,
            loadedDisplacement: SIMD2<Float>(0, -endDisplacement)
        )

        return FEMProblem2D(
            mesh: mesh,
            material: liangMaterial(),
            thickness: 1.0,
            integrationScheme: integrationScheme,
            prescribedDisplacements: prescribed,
            controls: benchmarkSolverControls(loadSteps: loadSteps)
        )
    }

    public static func rightEdgeDensity(
        problem: FEMProblem2D,
        densities: [Float]
    ) throws -> Float {
        let prepared = try problem.mesh.prepare()
        if prepared.elements.isEmpty {
            return 0
        }

        let bounds = meshBounds(problem.mesh)
        let xThreshold = bounds.min.x + 0.8 * (bounds.max.x - bounds.min.x)

        var weightedDensity: Float = 0
        var weightSum: Float = 0

        for (elementIndex, element) in prepared.elements.enumerated() {
            guard elementIndex < densities.count else {
                continue
            }
            let nodeIDs = element.nodeIDs
            let centroid = 0.25 * (
                prepared.nodes[Int(nodeIDs[0])] +
                prepared.nodes[Int(nodeIDs[1])] +
                prepared.nodes[Int(nodeIDs[2])] +
                prepared.nodes[Int(nodeIDs[3])]
            )

            if centroid.x >= xThreshold {
                let weight = max(1e-8, element.area)
                weightedDensity += weight * densities[elementIndex]
                weightSum += weight
            }
        }

        if weightSum <= 1e-8 {
            return densities.reduce(0, +) / Float(max(1, densities.count))
        }
        return weightedDensity / weightSum
    }
}

private func amirMaterial() -> MaterialParameters {
    MaterialParameters(
        youngsModulus: 85_000,
        poissonRatio: 0.30,
        yieldStress: 250,
        hardeningModulus: 540,
        damageOnset: 0.010,
        damageSlope: 0.20,
        damageCap: 0.96
    )
}

private func liangMaterial() -> MaterialParameters {
    MaterialParameters(
        youngsModulus: 72_000,
        poissonRatio: 0.30,
        yieldStress: 360,
        hardeningModulus: 420,
        damageOnset: 0.018,
        damageSlope: 0.30,
        damageCap: 0.95
    )
}

private func benchmarkSolverControls(loadSteps: Int) -> SolverControls {
    SolverControls(
        loadSteps: max(2, loadSteps),
        maxNewtonIterations: 40,
        residualTolerance: 2e-3,
        finiteDifferenceStep: 2e-4,
        lineSearchFloor: 1.0 / 128.0,
        linearSolverTolerance: 2e-5,
        linearSolverMaxIterations: 1_000
    )
}

private func filteredQuadMesh(
    mesh: Mesh2D,
    keepElement: (SIMD2<Float>) -> Bool
) -> Mesh2D {
    let quads = mesh.linearizedQuads(elements: mesh.elements)
    var kept: [SIMD4<UInt32>] = []
    kept.reserveCapacity(quads.count)

    for quad in quads {
        let centroid = 0.25 * (
            mesh.nodes[Int(quad[0])] +
            mesh.nodes[Int(quad[1])] +
            mesh.nodes[Int(quad[2])] +
            mesh.nodes[Int(quad[3])]
        )
        if keepElement(centroid) {
            kept.append(quad)
        }
    }

    var usedNodeIDs = Set<UInt32>()
    for quad in kept {
        usedNodeIDs.insert(quad[0])
        usedNodeIDs.insert(quad[1])
        usedNodeIDs.insert(quad[2])
        usedNodeIDs.insert(quad[3])
    }

    let sortedNodeIDs = usedNodeIDs.sorted()
    var remap: [UInt32: UInt32] = [:]
    remap.reserveCapacity(sortedNodeIDs.count)
    var nodes: [SIMD2<Float>] = []
    nodes.reserveCapacity(sortedNodeIDs.count)

    for (newID, oldID) in sortedNodeIDs.enumerated() {
        remap[oldID] = UInt32(newID)
        nodes.append(mesh.nodes[Int(oldID)])
    }

    let elements = kept.compactMap { quad -> Element2D? in
        guard
            let n0 = remap[quad[0]],
            let n1 = remap[quad[1]],
            let n2 = remap[quad[2]],
            let n3 = remap[quad[3]]
        else {
            return nil
        }
        return .quad4(SIMD4<UInt32>(n0, n1, n2, n3))
    }

    return Mesh2D(nodes: nodes, elements: elements)
}

private func makePrescribedDisplacements(
    mesh: Mesh2D,
    fixedNodePredicate: (SIMD2<Float>) -> Bool,
    loadedNodes: [Int],
    loadedDisplacement: SIMD2<Float>
) -> [PrescribedDisplacement2D] {
    var prescribed: [PrescribedDisplacement2D] = []
    prescribed.reserveCapacity(mesh.nodes.count * 2)

    var loadedSet = Set(loadedNodes)
    for (nodeID, point) in mesh.nodes.enumerated() {
        if fixedNodePredicate(point) {
            prescribed.append(PrescribedDisplacement2D(node: nodeID, component: 0, value: 0))
            prescribed.append(PrescribedDisplacement2D(node: nodeID, component: 1, value: 0))
            loadedSet.remove(nodeID)
        }
    }

    for nodeID in loadedSet.sorted() {
        if abs(loadedDisplacement.x) > 0 {
            prescribed.append(PrescribedDisplacement2D(node: nodeID, component: 0, value: loadedDisplacement.x))
        }
        if abs(loadedDisplacement.y) > 0 {
            prescribed.append(PrescribedDisplacement2D(node: nodeID, component: 1, value: loadedDisplacement.y))
        }
    }

    return prescribed
}

private func boundaryNodes(
    in mesh: Mesh2D,
    where predicate: (SIMD2<Float>) -> Bool,
    fallbackTarget: SIMD2<Float>,
    fallbackCount: Int
) -> [Int] {
    var nodes: [Int] = []
    nodes.reserveCapacity(16)

    for (index, point) in mesh.nodes.enumerated() where predicate(point) {
        nodes.append(index)
    }
    if !nodes.isEmpty {
        return nodes
    }

    let sortedByDistance = mesh.nodes.enumerated().sorted { lhs, rhs in
        let dL = simd_length_squared(lhs.element - fallbackTarget)
        let dR = simd_length_squared(rhs.element - fallbackTarget)
        return dL < dR
    }

    let count = max(1, min(max(1, fallbackCount), sortedByDistance.count))
    return Array(sortedByDistance.prefix(count).map(\.offset))
}

private func meshBounds(_ mesh: Mesh2D) -> (min: SIMD2<Float>, max: SIMD2<Float>) {
    var minPoint = SIMD2<Float>(Float.greatestFiniteMagnitude, Float.greatestFiniteMagnitude)
    var maxPoint = SIMD2<Float>(-Float.greatestFiniteMagnitude, -Float.greatestFiniteMagnitude)

    for point in mesh.nodes {
        minPoint = simd.min(minPoint, point)
        maxPoint = simd.max(maxPoint, point)
    }
    return (min: minPoint, max: maxPoint)
}

private func coordinateTolerance(bounds: (min: SIMD2<Float>, max: SIMD2<Float>)) -> Float {
    let span = max(bounds.max.x - bounds.min.x, bounds.max.y - bounds.min.y)
    return max(1e-6, 1e-5 * span)
}
