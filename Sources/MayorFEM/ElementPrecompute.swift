import Foundation
import simd

struct ElementGeometry {
    var nodeIDs: SIMD4<UInt32>
    var gradN0: SIMD4<Float>
    var gradN1: SIMD4<Float>
    var gradN2: SIMD4<Float>
    var gradN3: SIMD4<Float>
    var volume: Float

    func gradient(localIndex: Int) -> SIMD3<Float> {
        switch localIndex {
        case 0:
            return SIMD3<Float>(gradN0.x, gradN0.y, gradN0.z)
        case 1:
            return SIMD3<Float>(gradN1.x, gradN1.y, gradN1.z)
        case 2:
            return SIMD3<Float>(gradN2.x, gradN2.y, gradN2.z)
        case 3:
            return SIMD3<Float>(gradN3.x, gradN3.y, gradN3.z)
        default:
            preconditionFailure("Local node index must be in 0...3")
        }
    }
}

struct PreparedMesh {
    var nodes: [SIMD3<Float>]
    var elements: [ElementGeometry]
}

extension Mesh {
    func prepare() throws -> PreparedMesh {
        var preparedElements: [ElementGeometry] = []
        preparedElements.reserveCapacity(elements.count)

        for (elementIndex, connectivity) in elements.enumerated() {
            let ids = [Int(connectivity[0]), Int(connectivity[1]), Int(connectivity[2]), Int(connectivity[3])]
            for id in ids where id < 0 || id >= nodes.count {
                throw FEMError.invalidMesh("Element \(elementIndex) references out-of-range node \(id).")
            }

            var oriented = connectivity
            var elementGeometry = try buildElementGeometry(from: oriented)
            if elementGeometry.volume < 0 {
                oriented[1] = connectivity[2]
                oriented[2] = connectivity[1]
                elementGeometry = try buildElementGeometry(from: oriented)
            }

            if elementGeometry.volume <= 0 {
                throw FEMError.invalidMesh("Element \(elementIndex) has non-positive reference volume.")
            }

            preparedElements.append(elementGeometry)
        }

        return PreparedMesh(nodes: nodes, elements: preparedElements)
    }

    private func buildElementGeometry(from connectivity: SIMD4<UInt32>) throws -> ElementGeometry {
        let x1 = nodes[Int(connectivity[0])]
        let x2 = nodes[Int(connectivity[1])]
        let x3 = nodes[Int(connectivity[2])]
        let x4 = nodes[Int(connectivity[3])]

        let jacobian = simd_float3x3(columns: (x2 - x1, x3 - x1, x4 - x1))
        let detJ = simd_determinant(jacobian)

        if abs(detJ) < 1e-10 {
            throw FEMError.invalidMesh("Degenerate tetrahedron detected.")
        }

        let invJT = simd_transpose(simd_inverse(jacobian))
        let gradN1 = invJT * SIMD3<Float>(1, 0, 0)
        let gradN2 = invJT * SIMD3<Float>(0, 1, 0)
        let gradN3 = invJT * SIMD3<Float>(0, 0, 1)
        let gradN0 = -(gradN1 + gradN2 + gradN3)

        return ElementGeometry(
            nodeIDs: connectivity,
            gradN0: SIMD4<Float>(gradN0.x, gradN0.y, gradN0.z, 0),
            gradN1: SIMD4<Float>(gradN1.x, gradN1.y, gradN1.z, 0),
            gradN2: SIMD4<Float>(gradN2.x, gradN2.y, gradN2.z, 0),
            gradN3: SIMD4<Float>(gradN3.x, gradN3.y, gradN3.z, 0),
            volume: detJ / 6.0
        )
    }
}
