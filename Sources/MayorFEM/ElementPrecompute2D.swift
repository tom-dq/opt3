import Foundation
import simd

struct Quad4Geometry2D {
    var nodeIDs: SIMD4<UInt32>
    var area: Float

    func nodeID(_ index: Int) -> UInt32 {
        nodeIDs[index]
    }
}

struct PreparedMesh2D {
    var nodes: [SIMD2<Float>]
    var elements: [Quad4Geometry2D]
}

extension Mesh2D {
    func prepare() throws -> PreparedMesh2D {
        let quads = linearizedQuads(elements: elements)
        var prepared: [Quad4Geometry2D] = []
        prepared.reserveCapacity(quads.count)

        for (elementIndex, quad) in quads.enumerated() {
            let ids = [Int(quad[0]), Int(quad[1]), Int(quad[2]), Int(quad[3])]
            for id in ids where id < 0 || id >= nodes.count {
                throw FEMError.invalidMesh("2D element \(elementIndex) references out-of-range node \(id).")
            }

            var oriented = quad
            var area = polygonArea4(
                nodes[Int(oriented[0])],
                nodes[Int(oriented[1])],
                nodes[Int(oriented[2])],
                nodes[Int(oriented[3])]
            )

            if abs(area) < 1e-10 {
                throw FEMError.invalidMesh("2D quad element \(elementIndex) is degenerate.")
            }

            if area < 0 {
                oriented = SIMD4<UInt32>(quad[0], quad[3], quad[2], quad[1])
                area = -area
            }

            prepared.append(Quad4Geometry2D(nodeIDs: oriented, area: area))
        }

        return PreparedMesh2D(nodes: nodes, elements: prepared)
    }
}

@inline(__always)
private func polygonArea4(_ p0: SIMD2<Float>, _ p1: SIMD2<Float>, _ p2: SIMD2<Float>, _ p3: SIMD2<Float>) -> Float {
    0.5 * (
        p0.x * p1.y - p1.x * p0.y +
        p1.x * p2.y - p2.x * p1.y +
        p2.x * p3.y - p3.x * p2.y +
        p3.x * p0.y - p0.x * p3.y
    )
}
