import Foundation
import simd

public enum FEMVisualization {
    @discardableResult
    public static func writeVTKSeries(
        problem: FEMProblem,
        result: SolveResult,
        outputDirectory: String,
        deformationScale: Float = 5.0
    ) throws -> [String] {
        let preparedMesh = try problem.mesh.prepare()
        let snapshots = snapshotsForExport(result: result)
        let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)

        try FileManager.default.createDirectory(
            at: directoryURL,
            withIntermediateDirectories: true
        )

        var writtenFiles: [String] = []
        writtenFiles.reserveCapacity(snapshots.count)

        for snapshot in snapshots {
            let fileName = String(format: "step_%04d.vtk", snapshot.step)
            let fileURL = directoryURL.appendingPathComponent(fileName)

            let content = vtkContent(
                problem: problem,
                preparedMesh: preparedMesh,
                snapshot: snapshot,
                deformationScale: deformationScale
            )
            try content.write(to: fileURL, atomically: true, encoding: .utf8)
            writtenFiles.append(fileURL.path)
        }

        let pvdURL = directoryURL.appendingPathComponent("series.pvd")
        let pvdContent = pvdContent(for: snapshots)
        try pvdContent.write(to: pvdURL, atomically: true, encoding: .utf8)

        return writtenFiles
    }

    private static func snapshotsForExport(result: SolveResult) -> [StepSnapshot] {
        if !result.stepSnapshots.isEmpty {
            return result.stepSnapshots
        }

        return [
            StepSnapshot(
                step: 1,
                loadFactor: 1.0,
                displacements: result.displacements,
                elementStates: result.elementStates
            )
        ]
    }

    private static func vtkContent(
        problem: FEMProblem,
        preparedMesh: PreparedMesh,
        snapshot: StepSnapshot,
        deformationScale: Float
    ) -> String {
        let nodes = preparedMesh.nodes
        let elements = preparedMesh.elements
        let nodeCount = nodes.count
        let elementCount = elements.count

        let prescribed = prescribedNodeData(problem: problem, loadFactor: snapshot.loadFactor)
        let vonMises = elementVonMises(
            preparedMesh: preparedMesh,
            material: problem.material,
            displacements: snapshot.displacements,
            states: snapshot.elementStates
        )

        var lines: [String] = []
        lines.reserveCapacity(nodeCount * 8 + elementCount * 10 + 64)

        lines.append("# vtk DataFile Version 3.0")
        lines.append("MayorFEM step \(snapshot.step)")
        lines.append("ASCII")
        lines.append("DATASET UNSTRUCTURED_GRID")
        lines.append("POINTS \(nodeCount) float")

        for index in 0..<nodeCount {
            let p = nodes[index] + deformationScale * snapshot.displacements[index]
            lines.append("\(p.x) \(p.y) \(p.z)")
        }

        lines.append("CELLS \(elementCount) \(elementCount * 5)")
        for element in elements {
            lines.append("4 \(Int(element.nodeIDs[0])) \(Int(element.nodeIDs[1])) \(Int(element.nodeIDs[2])) \(Int(element.nodeIDs[3]))")
        }

        lines.append("CELL_TYPES \(elementCount)")
        for _ in 0..<elementCount {
            lines.append("10")
        }

        lines.append("POINT_DATA \(nodeCount)")
        lines.append("VECTORS displacement float")
        for displacement in snapshot.displacements {
            lines.append("\(displacement.x) \(displacement.y) \(displacement.z)")
        }

        lines.append("VECTORS prescribed_displacement float")
        for value in prescribed.vectors {
            lines.append("\(value.x) \(value.y) \(value.z)")
        }

        lines.append("SCALARS enforced_dof_count int 1")
        lines.append("LOOKUP_TABLE default")
        for count in prescribed.enforcedDOFCount {
            lines.append("\(count)")
        }

        lines.append("SCALARS prescribed_node int 1")
        lines.append("LOOKUP_TABLE default")
        for isPrescribed in prescribed.prescribedNode {
            lines.append("\(isPrescribed)")
        }

        lines.append("CELL_DATA \(elementCount)")
        lines.append("SCALARS von_mises float 1")
        lines.append("LOOKUP_TABLE default")
        for stress in vonMises {
            lines.append("\(stress)")
        }

        lines.append("SCALARS eq_plastic_strain float 1")
        lines.append("LOOKUP_TABLE default")
        for state in snapshot.elementStates {
            lines.append("\(state.equivalentPlasticStrain)")
        }

        lines.append("SCALARS damage float 1")
        lines.append("LOOKUP_TABLE default")
        for state in snapshot.elementStates {
            lines.append("\(state.damage)")
        }

        lines.append("SCALARS load_factor float 1")
        lines.append("LOOKUP_TABLE default")
        for _ in 0..<elementCount {
            lines.append("\(snapshot.loadFactor)")
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func prescribedNodeData(
        problem: FEMProblem,
        loadFactor: Float
    ) -> (vectors: [SIMD3<Float>], enforcedDOFCount: [Int], prescribedNode: [Int]) {
        let nodeCount = problem.mesh.nodes.count
        var vectors = Array(repeating: SIMD3<Float>.zero, count: nodeCount)
        var masks = Array(repeating: UInt8(0), count: nodeCount)

        for bc in problem.prescribedDisplacements {
            let scaledValue = loadFactor * bc.value
            vectors[bc.node][bc.component] = scaledValue
            masks[bc.node] |= UInt8(1 << bc.component)
        }

        var enforcedDOFCount = Array(repeating: 0, count: nodeCount)
        var prescribedNode = Array(repeating: 0, count: nodeCount)

        for index in 0..<nodeCount {
            let mask = masks[index]
            let count = Int(mask & 1) + Int((mask >> 1) & 1) + Int((mask >> 2) & 1)
            enforcedDOFCount[index] = count
            prescribedNode[index] = count > 0 ? 1 : 0
        }

        return (vectors: vectors, enforcedDOFCount: enforcedDOFCount, prescribedNode: prescribedNode)
    }

    private static func elementVonMises(
        preparedMesh: PreparedMesh,
        material: MaterialParameters,
        displacements: [SIMD3<Float>],
        states: [ElementState]
    ) -> [Float] {
        var vonMises: [Float] = []
        vonMises.reserveCapacity(preparedMesh.elements.count)

        for elementIndex in preparedMesh.elements.indices {
            let geometry = preparedMesh.elements[elementIndex]
            let reference = gatherElementNodalVectors(geometry: geometry, values: preparedMesh.nodes)
            let localDisplacements = gatherElementNodalVectors(geometry: geometry, values: displacements)
            let response = computeLocalElementResponse(
                geometry: geometry,
                reference: reference,
                displacements: localDisplacements,
                previousState: states[elementIndex],
                material: material
            )
            vonMises.append(response.vonMisesStress)
        }

        return vonMises
    }

    private static func pvdContent(for snapshots: [StepSnapshot]) -> String {
        var lines: [String] = []
        lines.append("<?xml version=\"1.0\"?>")
        lines.append("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
        lines.append("  <Collection>")

        for snapshot in snapshots {
            let fileName = String(format: "step_%04d.vtk", snapshot.step)
            lines.append("    <DataSet timestep=\"\(snapshot.loadFactor)\" group=\"\" part=\"0\" file=\"\(fileName)\"/>")
        }

        lines.append("  </Collection>")
        lines.append("</VTKFile>")
        return lines.joined(separator: "\n") + "\n"
    }
}
