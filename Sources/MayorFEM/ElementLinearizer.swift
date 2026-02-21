import Foundation
import simd

struct ElementLinearization {
    var nodeIDs: SIMD4<UInt32>
    var localResidual: [Float]
    var localTangent: [Float]
    var updatedState: ElementState
}

final class CPUElementLinearizer {
    private let preparedMesh: PreparedMesh
    private let material: MaterialParameters

    init(preparedMesh: PreparedMesh, material: MaterialParameters) {
        self.preparedMesh = preparedMesh
        self.material = material
    }

    func linearize(
        elementIndex: Int,
        displacements: [SIMD3<Float>],
        previousState: ElementState,
        finiteDifferenceStep: Float
    ) -> ElementLinearization {
        let geometry = preparedMesh.elements[elementIndex]
        let reference = gatherElementNodalVectors(geometry: geometry, values: preparedMesh.nodes)
        let localDisplacements = gatherElementNodalVectors(geometry: geometry, values: displacements)

        let baseResponse = computeLocalElementResponse(
            geometry: geometry,
            reference: reference,
            displacements: localDisplacements,
            previousState: previousState,
            material: material
        )
        let baseResidual = flattenElementNodalVectors(baseResponse.forces)
        var localTangent = Array(repeating: Float(0), count: 12 * 12)

        for column in 0..<12 {
            let node = column / 3
            let component = column % 3

            var perturbedDisplacements = localDisplacements
            var vector = perturbedDisplacements[node]
            let perturbation = finiteDifferenceStep * max(1.0, abs(vector[component]))
            vector[component] += perturbation
            perturbedDisplacements[node] = vector

            let perturbedResponse = computeLocalElementResponse(
                geometry: geometry,
                reference: reference,
                displacements: perturbedDisplacements,
                previousState: previousState,
                material: material
            )
            let perturbedResidual = flattenElementNodalVectors(perturbedResponse.forces)

            for row in 0..<12 {
                localTangent[row * 12 + column] =
                    (perturbedResidual[row] - baseResidual[row]) / perturbation
            }
        }

        return ElementLinearization(
            nodeIDs: geometry.nodeIDs,
            localResidual: baseResidual,
            localTangent: localTangent,
            updatedState: baseResponse.updatedState
        )
    }
}
