import Foundation
import simd
import CoreGraphics
import ImageIO

private enum Visualization2DError: Error {
    case invalidImageSize
    case contextCreationFailed
    case imageCreationFailed
    case pngDestinationCreationFailed
    case pngWriteFailed
}

extension FEMVisualization {
    @discardableResult
    public static func writeVTKSeries(
        problem: FEMProblem2D,
        result: SolveResult2D,
        outputDirectory: String,
        deformationScale: Float = 5.0
    ) throws -> [String] {
        let preparedMesh = try problem.mesh.prepare()
        let snapshots = result.stepSnapshots
        let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)

        try FileManager.default.createDirectory(
            at: directoryURL,
            withIntermediateDirectories: true
        )

        var files: [String] = []
        files.reserveCapacity(snapshots.count)

        for snapshot in snapshots {
            let fileName = String(format: "step_%04d.vtk", snapshot.step)
            let fileURL = directoryURL.appendingPathComponent(fileName)
            let content = vtkContent2D(
                problem: problem,
                preparedMesh: preparedMesh,
                snapshot: snapshot,
                deformationScale: deformationScale
            )
            try content.write(to: fileURL, atomically: true, encoding: .utf8)
            files.append(fileURL.path)
        }

        let pvdURL = directoryURL.appendingPathComponent("series.pvd")
        try pvdContent2D(for: snapshots).write(to: pvdURL, atomically: true, encoding: .utf8)

        return files
    }

    @discardableResult
    public static func writePNGSeries(
        problem: FEMProblem2D,
        result: SolveResult2D,
        outputDirectory: String,
        deformationScale: Float = 6.0,
        imageWidth: Int = 1400,
        imageHeight: Int = 900
    ) throws -> [String] {
        guard imageWidth > 0, imageHeight > 0 else {
            throw Visualization2DError.invalidImageSize
        }

        let preparedMesh = try problem.mesh.prepare()
        let snapshots = result.stepSnapshots
        let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)

        try FileManager.default.createDirectory(
            at: directoryURL,
            withIntermediateDirectories: true
        )

        var stressMin = Float.greatestFiniteMagnitude
        var stressMax = -Float.greatestFiniteMagnitude
        for snapshot in snapshots {
            if let localMin = snapshot.elementVonMises.min(), let localMax = snapshot.elementVonMises.max() {
                stressMin = min(stressMin, localMin)
                stressMax = max(stressMax, localMax)
            }
        }
        if !stressMin.isFinite || !stressMax.isFinite || stressMax < stressMin {
            stressMin = 0
            stressMax = 1
        }

        let elementCount = preparedMesh.elements.count
        var densityMin = Float.greatestFiniteMagnitude
        var densityMax = -Float.greatestFiniteMagnitude
        for snapshot in snapshots {
            let densities = snapshot.elementDensities.count == elementCount
                ? snapshot.elementDensities
                : Array(repeating: Float(1), count: elementCount)
            if let localMin = densities.min(), let localMax = densities.max() {
                densityMin = min(densityMin, localMin)
                densityMax = max(densityMax, localMax)
            }
        }
        if !densityMin.isFinite || !densityMax.isFinite || densityMax < densityMin {
            densityMin = 0
            densityMax = 1
        }

        var files: [String] = []
        files.reserveCapacity(snapshots.count)

        for snapshot in snapshots {
            let fileName = String(format: "step_%04d.png", snapshot.step)
            let fileURL = directoryURL.appendingPathComponent(fileName)
            try writePNGFrame2D(
                problem: problem,
                preparedMesh: preparedMesh,
                snapshot: snapshot,
                stressRange: (stressMin, stressMax),
                densityRange: (densityMin, densityMax),
                deformationScale: deformationScale,
                imageWidth: imageWidth,
                imageHeight: imageHeight,
                outputURL: fileURL
            )
            files.append(fileURL.path)
        }

        return files
    }
}

private func vtkContent2D(
    problem: FEMProblem2D,
    preparedMesh: PreparedMesh2D,
    snapshot: StepSnapshot2D,
    deformationScale: Float
) -> String {
    let nodeCount = preparedMesh.nodes.count
    let elementCount = preparedMesh.elements.count

    let prescribed = prescribedNodeData2D(problem: problem, loadFactor: snapshot.loadFactor)

    var lines: [String] = []
    lines.append("# vtk DataFile Version 3.0")
    lines.append("MayorFEM 2D step \(snapshot.step)")
    lines.append("ASCII")
    lines.append("DATASET UNSTRUCTURED_GRID")
    lines.append("POINTS \(nodeCount) float")

    for index in 0..<nodeCount {
        let p = preparedMesh.nodes[index] + deformationScale * snapshot.displacements[index]
        lines.append("\(p.x) \(p.y) 0")
    }

    lines.append("CELLS \(elementCount) \(elementCount * 5)")
    for element in preparedMesh.elements {
        lines.append("4 \(Int(element.nodeIDs[0])) \(Int(element.nodeIDs[1])) \(Int(element.nodeIDs[2])) \(Int(element.nodeIDs[3]))")
    }

    lines.append("CELL_TYPES \(elementCount)")
    for _ in 0..<elementCount {
        lines.append("9")
    }

    lines.append("POINT_DATA \(nodeCount)")
    lines.append("VECTORS displacement float")
    for displacement in snapshot.displacements {
        lines.append("\(displacement.x) \(displacement.y) 0")
    }

    lines.append("VECTORS prescribed_displacement float")
    for value in prescribed.vectors {
        lines.append("\(value.x) \(value.y) 0")
    }

    lines.append("SCALARS enforced_dof_count int 1")
    lines.append("LOOKUP_TABLE default")
    for count in prescribed.enforcedDOFCount {
        lines.append("\(count)")
    }

    lines.append("SCALARS prescribed_node int 1")
    lines.append("LOOKUP_TABLE default")
    for flag in prescribed.prescribedNode {
        lines.append("\(flag)")
    }

    lines.append("CELL_DATA \(elementCount)")
    lines.append("SCALARS von_mises float 1")
    lines.append("LOOKUP_TABLE default")
    for stress in snapshot.elementVonMises {
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

    let densities = snapshot.elementDensities.count == elementCount
        ? snapshot.elementDensities
        : Array(repeating: Float(1), count: elementCount)
    lines.append("SCALARS density float 1")
    lines.append("LOOKUP_TABLE default")
    for density in densities {
        lines.append("\(density)")
    }

    return lines.joined(separator: "\n") + "\n"
}

private func pvdContent2D(for snapshots: [StepSnapshot2D]) -> String {
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

private func prescribedNodeData2D(
    problem: FEMProblem2D,
    loadFactor: Float
) -> (vectors: [SIMD2<Float>], enforcedDOFCount: [Int], prescribedNode: [Int]) {
    let nodeCount = problem.mesh.nodes.count
    var vectors = Array(repeating: SIMD2<Float>.zero, count: nodeCount)
    var masks = Array(repeating: UInt8(0), count: nodeCount)

    for bc in problem.prescribedDisplacements {
        let scaled = loadFactor * bc.value
        vectors[bc.node][bc.component] = scaled
        masks[bc.node] |= UInt8(1 << bc.component)
    }

    var count = Array(repeating: 0, count: nodeCount)
    var nodeFlag = Array(repeating: 0, count: nodeCount)

    for index in 0..<nodeCount {
        let mask = masks[index]
        let c = Int(mask & 1) + Int((mask >> 1) & 1)
        count[index] = c
        nodeFlag[index] = c > 0 ? 1 : 0
    }

    return (vectors: vectors, enforcedDOFCount: count, prescribedNode: nodeFlag)
}

private func writePNGFrame2D(
    problem: FEMProblem2D,
    preparedMesh: PreparedMesh2D,
    snapshot: StepSnapshot2D,
    stressRange: (Float, Float),
    densityRange: (Float, Float),
    deformationScale: Float,
    imageWidth: Int,
    imageHeight: Int,
    outputURL: URL
) throws {
    let points = deformedPoints2D(
        nodes: preparedMesh.nodes,
        displacements: snapshot.displacements,
        deformationScale: deformationScale
    )

    let imageWidthF = Float(imageWidth)
    let imageHeightF = Float(imageHeight)
    let outerMargin: Float = 22
    let panelGap: Float = 20
    let panelWidth = max(1, 0.5 * (imageWidthF - 2 * outerMargin - panelGap))
    let panelHeight = max(1, imageHeightF - 2 * outerMargin)

    let stressViewport = Viewport2D(
        x: outerMargin,
        y: outerMargin,
        width: panelWidth,
        height: panelHeight
    )
    let densityViewport = Viewport2D(
        x: outerMargin + panelWidth + panelGap,
        y: outerMargin,
        width: panelWidth,
        height: panelHeight
    )

    let stressFrame = frameTransform2D(points: points, imageHeight: imageHeight, viewport: stressViewport)
    let densityFrame = frameTransform2D(points: points, imageHeight: imageHeight, viewport: densityViewport)

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
        throw Visualization2DError.contextCreationFailed
    }

    context.setFillColor(CGColor(red: 0.97, green: 0.98, blue: 0.99, alpha: 1))
    context.fill(CGRect(x: 0, y: 0, width: imageWidth, height: imageHeight))

    context.setFillColor(CGColor(red: 1, green: 1, blue: 1, alpha: 0.98))
    context.fill(CGRect(
        x: CGFloat(stressViewport.x),
        y: CGFloat(stressViewport.y),
        width: CGFloat(stressViewport.width),
        height: CGFloat(stressViewport.height)
    ))
    context.fill(CGRect(
        x: CGFloat(densityViewport.x),
        y: CGFloat(densityViewport.y),
        width: CGFloat(densityViewport.width),
        height: CGFloat(densityViewport.height)
    ))

    let densities = snapshot.elementDensities.count == preparedMesh.elements.count
        ? snapshot.elementDensities
        : Array(repeating: Float(1), count: preparedMesh.elements.count)

    for (elementIndex, element) in preparedMesh.elements.enumerated() {
        let ids = element.nodeIDs
        let s0 = stressFrame.project(points[Int(ids[0])])
        let s1 = stressFrame.project(points[Int(ids[1])])
        let s2 = stressFrame.project(points[Int(ids[2])])
        let s3 = stressFrame.project(points[Int(ids[3])])
        let d0 = densityFrame.project(points[Int(ids[0])])
        let d1 = densityFrame.project(points[Int(ids[1])])
        let d2 = densityFrame.project(points[Int(ids[2])])
        let d3 = densityFrame.project(points[Int(ids[3])])

        let stress = snapshot.elementVonMises[elementIndex]
        let t = normalize2D(value: stress, minValue: stressRange.0, maxValue: stressRange.1)
        let color = heatMapColor2D(t)

        context.beginPath()
        context.move(to: s0)
        context.addLine(to: s1)
        context.addLine(to: s2)
        context.addLine(to: s3)
        context.closePath()
        context.setFillColor(CGColor(
            red: CGFloat(color.x),
            green: CGFloat(color.y),
            blue: CGFloat(color.z),
            alpha: 1
        ))
        context.fillPath()

        context.beginPath()
        context.move(to: s0)
        context.addLine(to: s1)
        context.addLine(to: s2)
        context.addLine(to: s3)
        context.closePath()
        context.setStrokeColor(CGColor(red: 0, green: 0, blue: 0, alpha: 0.25))
        context.setLineWidth(0.6)
        context.strokePath()

        let density = densities[elementIndex]
        let densityT = normalize2D(value: density, minValue: densityRange.0, maxValue: densityRange.1)
        let densityColor = densityMapColor2D(densityT)

        context.beginPath()
        context.move(to: d0)
        context.addLine(to: d1)
        context.addLine(to: d2)
        context.addLine(to: d3)
        context.closePath()
        context.setFillColor(CGColor(
            red: CGFloat(densityColor.x),
            green: CGFloat(densityColor.y),
            blue: CGFloat(densityColor.z),
            alpha: 1
        ))
        context.fillPath()

        context.beginPath()
        context.move(to: d0)
        context.addLine(to: d1)
        context.addLine(to: d2)
        context.addLine(to: d3)
        context.closePath()
        context.setStrokeColor(CGColor(red: 0, green: 0, blue: 0, alpha: 0.25))
        context.setLineWidth(0.6)
        context.strokePath()
    }

    drawPrescribedArrows2D(
        context: context,
        problem: problem,
        snapshot: snapshot,
        points: points,
        frame: stressFrame,
        deformationScale: deformationScale
    )
    drawPrescribedArrows2D(
        context: context,
        problem: problem,
        snapshot: snapshot,
        points: points,
        frame: densityFrame,
        deformationScale: deformationScale
    )

    context.setStrokeColor(CGColor(red: 0.15, green: 0.15, blue: 0.18, alpha: 0.45))
    context.setLineWidth(1.0)
    context.stroke(CGRect(
        x: CGFloat(stressViewport.x),
        y: CGFloat(stressViewport.y),
        width: CGFloat(stressViewport.width),
        height: CGFloat(stressViewport.height)
    ))
    context.stroke(CGRect(
        x: CGFloat(densityViewport.x),
        y: CGFloat(densityViewport.y),
        width: CGFloat(densityViewport.width),
        height: CGFloat(densityViewport.height)
    ))

    guard let image = context.makeImage() else {
        throw Visualization2DError.imageCreationFailed
    }

    guard let destination = CGImageDestinationCreateWithURL(
        outputURL as CFURL,
        "public.png" as CFString,
        1,
        nil
    ) else {
        throw Visualization2DError.pngDestinationCreationFailed
    }

    CGImageDestinationAddImage(destination, image, nil)
    guard CGImageDestinationFinalize(destination) else {
        throw Visualization2DError.pngWriteFailed
    }
}

private struct FrameTransform2D {
    var minX: Float
    var minY: Float
    var scale: Float
    var imageHeight: Float
    var viewport: Viewport2D
    var margin: Float

    func project(_ point: SIMD2<Float>) -> CGPoint {
        let x = viewport.x + margin + (point.x - minX) * scale
        let y = imageHeight - (viewport.y + margin + (point.y - minY) * scale)
        return CGPoint(x: CGFloat(x), y: CGFloat(y))
    }
}

private struct Viewport2D {
    var x: Float
    var y: Float
    var width: Float
    var height: Float
}

private func deformedPoints2D(
    nodes: [SIMD2<Float>],
    displacements: [SIMD2<Float>],
    deformationScale: Float
) -> [SIMD2<Float>] {
    var points = Array(repeating: SIMD2<Float>.zero, count: nodes.count)
    for index in nodes.indices {
        points[index] = nodes[index] + deformationScale * displacements[index]
    }
    return points
}

private func frameTransform2D(points: [SIMD2<Float>], imageHeight: Int, viewport: Viewport2D) -> FrameTransform2D {
    var minX = Float.greatestFiniteMagnitude
    var minY = Float.greatestFiniteMagnitude
    var maxX = -Float.greatestFiniteMagnitude
    var maxY = -Float.greatestFiniteMagnitude

    for point in points {
        minX = min(minX, point.x)
        minY = min(minY, point.y)
        maxX = max(maxX, point.x)
        maxY = max(maxY, point.y)
    }

    let margin: Float = 14
    let spanX = max(1e-6, maxX - minX)
    let spanY = max(1e-6, maxY - minY)

    let usableWidth = max(1, viewport.width - 2 * margin)
    let usableHeight = max(1, viewport.height - 2 * margin)
    let scale = min(usableWidth / spanX, usableHeight / spanY)

    return FrameTransform2D(
        minX: minX,
        minY: minY,
        scale: scale,
        imageHeight: Float(imageHeight),
        viewport: viewport,
        margin: margin
    )
}

private func drawPrescribedArrows2D(
    context: CGContext,
    problem: FEMProblem2D,
    snapshot: StepSnapshot2D,
    points: [SIMD2<Float>],
    frame: FrameTransform2D,
    deformationScale: Float
) {
    let prescribed = prescribedNodeData2D(problem: problem, loadFactor: snapshot.loadFactor)

    for nodeIndex in problem.mesh.nodes.indices {
        if prescribed.enforcedDOFCount[nodeIndex] == 0 {
            continue
        }

        let start = frame.project(points[nodeIndex])
        let v = prescribed.vectors[nodeIndex]

        context.setFillColor(CGColor(red: 0.85, green: 0.12, blue: 0.12, alpha: 0.9))
        let radius: CGFloat = 3.2
        context.fillEllipse(in: CGRect(x: start.x - radius, y: start.y - radius, width: 2 * radius, height: 2 * radius))

        let norm = simd_length(v)
        if norm < 1e-8 {
            continue
        }

        let endWorld = points[nodeIndex] + 1.6 * deformationScale * v
        let end = frame.project(endWorld)

        context.setStrokeColor(CGColor(red: 0.78, green: 0.12, blue: 0.12, alpha: 0.95))
        context.setLineWidth(1.5)
        context.beginPath()
        context.move(to: start)
        context.addLine(to: end)
        context.strokePath()
    }
}

private func normalize2D(value: Float, minValue: Float, maxValue: Float) -> Float {
    let span = maxValue - minValue
    if span < 1e-9 {
        return 0.5
    }
    return Swift.max(0, Swift.min(1, (value - minValue) / span))
}

private func heatMapColor2D(_ t: Float) -> SIMD3<Float> {
    let x = Swift.max(0, Swift.min(1, t))

    if x < 0.25 {
        let s = x / 0.25
        return SIMD3<Float>(0.1, 0.2 + 0.7 * s, 0.8 + 0.2 * s)
    } else if x < 0.5 {
        let s = (x - 0.25) / 0.25
        return SIMD3<Float>(0.1 + 0.3 * s, 0.9, 1.0 - 0.7 * s)
    } else if x < 0.75 {
        let s = (x - 0.5) / 0.25
        return SIMD3<Float>(0.4 + 0.55 * s, 0.9 - 0.2 * s, 0.3 - 0.2 * s)
    } else {
        let s = (x - 0.75) / 0.25
        return SIMD3<Float>(0.95, 0.7 - 0.45 * s, 0.1)
    }
}

private func densityMapColor2D(_ t: Float) -> SIMD3<Float> {
    let x = Swift.max(0, Swift.min(1, t))
    if x < 0.5 {
        let s = x / 0.5
        return SIMD3<Float>(0.08 + 0.32 * s, 0.12 + 0.58 * s, 0.38 + 0.47 * s)
    } else {
        let s = (x - 0.5) / 0.5
        return SIMD3<Float>(0.40 + 0.55 * s, 0.70 - 0.18 * s, 0.85 - 0.62 * s)
    }
}
