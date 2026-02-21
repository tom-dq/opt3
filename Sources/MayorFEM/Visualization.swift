import Foundation
import simd
import CoreGraphics
import ImageIO

private enum VisualizationError: Error {
    case invalidImageSize
    case contextCreationFailed
    case imageCreationFailed
    case pngDestinationCreationFailed
    case pngWriteFailed
}

private struct SurfaceFace {
    var nodeIDs: SIMD3<Int>
    var elementIndex: Int
}

private struct ProjectedTriangle {
    var p0: CGPoint
    var p1: CGPoint
    var p2: CGPoint
    var depth: Float
    var stress: Float
    var shade: Float
}

private struct ProjectionTransform {
    var center: SIMD3<Float>
    var rotation: simd_float3x3
    var minX: Float
    var minY: Float
    var scale: Float
    var imageHeight: Float
    var margin: Float

    func cameraPoint(_ world: SIMD3<Float>) -> SIMD3<Float> {
        rotation * (world - center)
    }

    func screenPoint(_ world: SIMD3<Float>) -> CGPoint {
        let p = cameraPoint(world)
        let x = margin + (p.x - minX) * scale
        let y = imageHeight - (margin + (p.y - minY) * scale)
        return CGPoint(x: CGFloat(x), y: CGFloat(y))
    }
}

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

    @discardableResult
    public static func writePNGSeries(
        problem: FEMProblem,
        result: SolveResult,
        outputDirectory: String,
        deformationScale: Float = 5.0,
        imageWidth: Int = 1280,
        imageHeight: Int = 900
    ) throws -> [String] {
        guard imageWidth > 0, imageHeight > 0 else {
            throw VisualizationError.invalidImageSize
        }

        let preparedMesh = try problem.mesh.prepare()
        let snapshots = snapshotsForExport(result: result)
        let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)

        try FileManager.default.createDirectory(
            at: directoryURL,
            withIntermediateDirectories: true
        )

        let surfaceFaces = boundaryFaces(of: preparedMesh)

        var perSnapshotStress: [[Float]] = []
        perSnapshotStress.reserveCapacity(snapshots.count)

        var globalStressMin = Float.greatestFiniteMagnitude
        var globalStressMax = -Float.greatestFiniteMagnitude

        for snapshot in snapshots {
            let stress = elementVonMises(
                preparedMesh: preparedMesh,
                material: problem.material,
                displacements: snapshot.displacements,
                states: snapshot.elementStates
            )
            perSnapshotStress.append(stress)

            if let localMin = stress.min(), let localMax = stress.max() {
                globalStressMin = min(globalStressMin, localMin)
                globalStressMax = max(globalStressMax, localMax)
            }
        }

        if !globalStressMin.isFinite || !globalStressMax.isFinite || globalStressMax < globalStressMin {
            globalStressMin = 0
            globalStressMax = 1
        }

        var writtenFiles: [String] = []
        writtenFiles.reserveCapacity(snapshots.count)

        for (snapshotIndex, snapshot) in snapshots.enumerated() {
            let fileName = String(format: "step_%04d.png", snapshot.step)
            let fileURL = directoryURL.appendingPathComponent(fileName)

            try writeSnapshotPNG(
                problem: problem,
                preparedMesh: preparedMesh,
                surfaceFaces: surfaceFaces,
                snapshot: snapshot,
                stress: perSnapshotStress[snapshotIndex],
                stressRange: (globalStressMin, globalStressMax),
                deformationScale: deformationScale,
                imageWidth: imageWidth,
                imageHeight: imageHeight,
                outputURL: fileURL
            )
            writtenFiles.append(fileURL.path)
        }

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

    private static func writeSnapshotPNG(
        problem: FEMProblem,
        preparedMesh: PreparedMesh,
        surfaceFaces: [SurfaceFace],
        snapshot: StepSnapshot,
        stress: [Float],
        stressRange: (Float, Float),
        deformationScale: Float,
        imageWidth: Int,
        imageHeight: Int,
        outputURL: URL
    ) throws {
        let worldPoints = deformedPoints(
            nodes: preparedMesh.nodes,
            displacements: snapshot.displacements,
            deformationScale: deformationScale
        )

        let transform = projectionTransform(
            worldPoints: worldPoints,
            imageWidth: imageWidth,
            imageHeight: imageHeight
        )

        let colorSpace = CGColorSpaceCreateDeviceRGB()
        guard let context = CGContext(
            data: nil,
            width: imageWidth,
            height: imageHeight,
            bitsPerComponent: 8,
            bytesPerRow: imageWidth * 4,
            space: colorSpace,
            bitmapInfo: CGImageAlphaInfo.premultipliedLast.rawValue
        ) else {
            throw VisualizationError.contextCreationFailed
        }

        context.setFillColor(CGColor(red: 0.96, green: 0.97, blue: 0.98, alpha: 1))
        context.fill(CGRect(x: 0, y: 0, width: imageWidth, height: imageHeight))

        let triangles = projectedTriangles(
            surfaceFaces: surfaceFaces,
            worldPoints: worldPoints,
            stress: stress,
            transform: transform
        )

        for triangle in triangles.sorted(by: { $0.depth < $1.depth }) {
            let t = normalizedStress(
                value: triangle.stress,
                minValue: stressRange.0,
                maxValue: stressRange.1
            )
            let base = heatMapColor(t)
            let shaded = SIMD4<Float>(
                min(1, base.x * triangle.shade),
                min(1, base.y * triangle.shade),
                min(1, base.z * triangle.shade),
                1
            )

            context.beginPath()
            context.move(to: triangle.p0)
            context.addLine(to: triangle.p1)
            context.addLine(to: triangle.p2)
            context.closePath()
            context.setFillColor(CGColor(
                red: CGFloat(shaded.x),
                green: CGFloat(shaded.y),
                blue: CGFloat(shaded.z),
                alpha: CGFloat(shaded.w)
            ))
            context.fillPath()

            context.beginPath()
            context.move(to: triangle.p0)
            context.addLine(to: triangle.p1)
            context.addLine(to: triangle.p2)
            context.closePath()
            context.setStrokeColor(CGColor(red: 0, green: 0, blue: 0, alpha: 0.15))
            context.setLineWidth(0.7)
            context.strokePath()
        }

        drawEnforcedDisplacements(
            context: context,
            problem: problem,
            snapshot: snapshot,
            worldPoints: worldPoints,
            transform: transform,
            deformationScale: deformationScale
        )

        guard let image = context.makeImage() else {
            throw VisualizationError.imageCreationFailed
        }

        guard let destination = CGImageDestinationCreateWithURL(
            outputURL as CFURL,
            "public.png" as CFString,
            1,
            nil
        ) else {
            throw VisualizationError.pngDestinationCreationFailed
        }

        CGImageDestinationAddImage(destination, image, nil)
        guard CGImageDestinationFinalize(destination) else {
            throw VisualizationError.pngWriteFailed
        }
    }

    private static func drawEnforcedDisplacements(
        context: CGContext,
        problem: FEMProblem,
        snapshot: StepSnapshot,
        worldPoints: [SIMD3<Float>],
        transform: ProjectionTransform,
        deformationScale: Float
    ) {
        let prescribed = prescribedNodeData(problem: problem, loadFactor: snapshot.loadFactor)

        for nodeIndex in problem.mesh.nodes.indices {
            if prescribed.enforcedDOFCount[nodeIndex] == 0 {
                continue
            }

            let start = transform.screenPoint(worldPoints[nodeIndex])
            let vector = prescribed.vectors[nodeIndex]

            context.setFillColor(CGColor(red: 0.82, green: 0.1, blue: 0.1, alpha: 0.95))
            let radius: CGFloat = 3.5
            context.fillEllipse(in: CGRect(
                x: start.x - radius,
                y: start.y - radius,
                width: 2 * radius,
                height: 2 * radius
            ))

            let vectorNorm = simd_length(vector)
            if vectorNorm < 1e-8 {
                continue
            }

            let endWorld = worldPoints[nodeIndex] + (1.5 * deformationScale) * vector
            let end = transform.screenPoint(endWorld)
            drawArrow(context: context, start: start, end: end)
        }
    }

    private static func drawArrow(context: CGContext, start: CGPoint, end: CGPoint) {
        let dx = end.x - start.x
        let dy = end.y - start.y
        let length = hypot(dx, dy)
        if length < 2 {
            return
        }

        let ux = dx / length
        let uy = dy / length

        context.setStrokeColor(CGColor(red: 0.78, green: 0.12, blue: 0.12, alpha: 0.95))
        context.setLineWidth(1.6)
        context.beginPath()
        context.move(to: start)
        context.addLine(to: end)
        context.strokePath()

        let headLength: CGFloat = min(10, max(6, 0.16 * length))
        let headWidth: CGFloat = headLength * 0.55

        let left = CGPoint(
            x: end.x - headLength * ux + headWidth * (-uy),
            y: end.y - headLength * uy + headWidth * ux
        )
        let right = CGPoint(
            x: end.x - headLength * ux - headWidth * (-uy),
            y: end.y - headLength * uy - headWidth * ux
        )

        context.beginPath()
        context.move(to: end)
        context.addLine(to: left)
        context.addLine(to: right)
        context.closePath()
        context.setFillColor(CGColor(red: 0.78, green: 0.12, blue: 0.12, alpha: 0.95))
        context.fillPath()
    }

    private static func projectedTriangles(
        surfaceFaces: [SurfaceFace],
        worldPoints: [SIMD3<Float>],
        stress: [Float],
        transform: ProjectionTransform
    ) -> [ProjectedTriangle] {
        let light = simd_normalize(SIMD3<Float>(0.35, 0.4, 1.0))

        var triangles: [ProjectedTriangle] = []
        triangles.reserveCapacity(surfaceFaces.count)

        for face in surfaceFaces {
            let p0World = worldPoints[face.nodeIDs.x]
            let p1World = worldPoints[face.nodeIDs.y]
            let p2World = worldPoints[face.nodeIDs.z]

            let c0 = transform.cameraPoint(p0World)
            let c1 = transform.cameraPoint(p1World)
            let c2 = transform.cameraPoint(p2World)

            let normal = simd_cross(c1 - c0, c2 - c0)
            let normalNorm = simd_length(normal)
            if normalNorm < 1e-8 {
                continue
            }

            let unitNormal = normal / normalNorm
            let shade = 0.35 + 0.65 * abs(simd_dot(unitNormal, light))

            triangles.append(
                ProjectedTriangle(
                    p0: transform.screenPoint(p0World),
                    p1: transform.screenPoint(p1World),
                    p2: transform.screenPoint(p2World),
                    depth: (c0.z + c1.z + c2.z) / 3,
                    stress: stress[face.elementIndex],
                    shade: shade
                )
            )
        }

        return triangles
    }

    private static func projectionTransform(
        worldPoints: [SIMD3<Float>],
        imageWidth: Int,
        imageHeight: Int
    ) -> ProjectionTransform {
        let center = worldPoints.reduce(SIMD3<Float>.zero, +) / Float(max(1, worldPoints.count))

        let yaw = Float.pi / 4.5
        let pitch = Float.pi / 7.0

        let rotationZ = simd_float3x3(
            SIMD3<Float>(cos(yaw), -sin(yaw), 0),
            SIMD3<Float>(sin(yaw), cos(yaw), 0),
            SIMD3<Float>(0, 0, 1)
        )
        let rotationX = simd_float3x3(
            SIMD3<Float>(1, 0, 0),
            SIMD3<Float>(0, cos(pitch), -sin(pitch)),
            SIMD3<Float>(0, sin(pitch), cos(pitch))
        )

        let rotation = rotationX * rotationZ
        let cameraPoints = worldPoints.map { rotation * ($0 - center) }

        var minX = Float.greatestFiniteMagnitude
        var maxX = -Float.greatestFiniteMagnitude
        var minY = Float.greatestFiniteMagnitude
        var maxY = -Float.greatestFiniteMagnitude

        for p in cameraPoints {
            minX = min(minX, p.x)
            maxX = max(maxX, p.x)
            minY = min(minY, p.y)
            maxY = max(maxY, p.y)
        }

        let margin: Float = 36
        let usableWidth = max(1, Float(imageWidth) - 2 * margin)
        let usableHeight = max(1, Float(imageHeight) - 2 * margin)

        let spanX = max(1e-6, maxX - minX)
        let spanY = max(1e-6, maxY - minY)
        let scale = min(usableWidth / spanX, usableHeight / spanY)

        return ProjectionTransform(
            center: center,
            rotation: rotation,
            minX: minX,
            minY: minY,
            scale: scale,
            imageHeight: Float(imageHeight),
            margin: margin
        )
    }

    private static func deformedPoints(
        nodes: [SIMD3<Float>],
        displacements: [SIMD3<Float>],
        deformationScale: Float
    ) -> [SIMD3<Float>] {
        var points = Array(repeating: SIMD3<Float>.zero, count: nodes.count)
        for index in nodes.indices {
            points[index] = nodes[index] + deformationScale * displacements[index]
        }
        return points
    }

    private static func boundaryFaces(of preparedMesh: PreparedMesh) -> [SurfaceFace] {
        struct FaceRecord {
            var count: Int
            var nodeIDs: SIMD3<Int>
            var elementIndex: Int
        }

        var records: [String: FaceRecord] = [:]

        for (elementIndex, element) in preparedMesh.elements.enumerated() {
            let n0 = Int(element.nodeIDs[0])
            let n1 = Int(element.nodeIDs[1])
            let n2 = Int(element.nodeIDs[2])
            let n3 = Int(element.nodeIDs[3])

            let faces = [
                SIMD3<Int>(n0, n1, n2),
                SIMD3<Int>(n0, n1, n3),
                SIMD3<Int>(n0, n2, n3),
                SIMD3<Int>(n1, n2, n3),
            ]

            for face in faces {
                let key = faceKey(face)
                if var record = records[key] {
                    record.count += 1
                    records[key] = record
                } else {
                    records[key] = FaceRecord(count: 1, nodeIDs: face, elementIndex: elementIndex)
                }
            }
        }

        var boundary: [SurfaceFace] = []
        boundary.reserveCapacity(records.count)

        for record in records.values where record.count == 1 {
            boundary.append(SurfaceFace(nodeIDs: record.nodeIDs, elementIndex: record.elementIndex))
        }

        return boundary
    }

    private static func faceKey(_ face: SIMD3<Int>) -> String {
        let sorted = [face.x, face.y, face.z].sorted()
        return "\(sorted[0])_\(sorted[1])_\(sorted[2])"
    }

    private static func normalizedStress(value: Float, minValue: Float, maxValue: Float) -> Float {
        let span = maxValue - minValue
        if span < 1e-9 {
            return 0.5
        }
        return Swift.max(0, Swift.min(1, (value - minValue) / span))
    }

    private static func heatMapColor(_ t: Float) -> SIMD4<Float> {
        let x = max(0, min(1, t))

        let r: Float
        let g: Float
        let b: Float

        if x < 0.25 {
            let s = x / 0.25
            r = 0.1
            g = 0.2 + 0.7 * s
            b = 0.8 + 0.2 * s
        } else if x < 0.5 {
            let s = (x - 0.25) / 0.25
            r = 0.1 + 0.3 * s
            g = 0.9
            b = 1.0 - 0.7 * s
        } else if x < 0.75 {
            let s = (x - 0.5) / 0.25
            r = 0.4 + 0.55 * s
            g = 0.9 - 0.2 * s
            b = 0.3 - 0.2 * s
        } else {
            let s = (x - 0.75) / 0.25
            r = 0.95
            g = 0.7 - 0.45 * s
            b = 0.1
        }

        return SIMD4<Float>(r, g, b, 1)
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
