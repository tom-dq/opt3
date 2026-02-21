import Foundation
import simd

struct Tri3Geometry2D {
    var nodeIDs: SIMD3<UInt32>
    var gradN0: SIMD2<Float>
    var gradN1: SIMD2<Float>
    var gradN2: SIMD2<Float>
    var area: Float

    func gradient(localIndex: Int) -> SIMD2<Float> {
        switch localIndex {
        case 0: return gradN0
        case 1: return gradN1
        case 2: return gradN2
        default:
            preconditionFailure("Local node index must be in 0...2")
        }
    }
}

struct PreparedMesh2D {
    var nodes: [SIMD2<Float>]
    var elements: [Tri3Geometry2D]
}

extension Mesh2D {
    func prepare() throws -> PreparedMesh2D {
        let triangles = linearizedTriangles(elements: elements)
        var prepared: [Tri3Geometry2D] = []
        prepared.reserveCapacity(triangles.count)

        for (elementIndex, triangle) in triangles.enumerated() {
            let ids = [Int(triangle.x), Int(triangle.y), Int(triangle.z)]
            for id in ids where id < 0 || id >= nodes.count {
                throw FEMError.invalidMesh("2D element \(elementIndex) references out-of-range node \(id).")
            }

            let x0 = nodes[ids[0]]
            let x1 = nodes[ids[1]]
            let x2 = nodes[ids[2]]

            let dx1 = x1.x - x0.x
            let dy1 = x1.y - x0.y
            let dx2 = x2.x - x0.x
            let dy2 = x2.y - x0.y

            let twoArea = dx1 * dy2 - dy1 * dx2
            if abs(twoArea) < 1e-10 {
                throw FEMError.invalidMesh("2D element \(elementIndex) is degenerate.")
            }

            var oriented = triangle
            var area = 0.5 * twoArea
            if area < 0 {
                oriented = SIMD3<UInt32>(triangle.x, triangle.z, triangle.y)
                area = -area
            }

            let a = Int(oriented.x)
            let b = Int(oriented.y)
            let c = Int(oriented.z)

            let xa = nodes[a]
            let xb = nodes[b]
            let xc = nodes[c]

            let denom = (xb.x - xa.x) * (xc.y - xa.y) - (xb.y - xa.y) * (xc.x - xa.x)
            if abs(denom) < 1e-10 {
                throw FEMError.invalidMesh("2D element \(elementIndex) has invalid orientation/area.")
            }

            let invDenom = 1.0 / denom

            let gradN0 = SIMD2<Float>(
                (xb.y - xc.y) * invDenom,
                (xc.x - xb.x) * invDenom
            )
            let gradN1 = SIMD2<Float>(
                (xc.y - xa.y) * invDenom,
                (xa.x - xc.x) * invDenom
            )
            let gradN2 = SIMD2<Float>(
                (xa.y - xb.y) * invDenom,
                (xb.x - xa.x) * invDenom
            )

            prepared.append(
                Tri3Geometry2D(
                    nodeIDs: oriented,
                    gradN0: gradN0,
                    gradN1: gradN1,
                    gradN2: gradN2,
                    area: area
                )
            )
        }

        return PreparedMesh2D(nodes: nodes, elements: prepared)
    }
}
