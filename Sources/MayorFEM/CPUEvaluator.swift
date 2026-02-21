import Foundation
import simd

final class CPUEvaluator: ElementEvaluator {
    private let preparedMesh: PreparedMesh
    private let material: MaterialParameters

    init(preparedMesh: PreparedMesh, material: MaterialParameters) {
        self.preparedMesh = preparedMesh
        self.material = material
    }

    func evaluate(displacements: [SIMD3<Float>], previousStates: [ElementState]) throws -> ElementEvaluation {
        if displacements.count != preparedMesh.nodes.count {
            throw FEMError.invalidMesh("Displacement vector length does not match node count.")
        }
        if previousStates.count != preparedMesh.elements.count {
            throw FEMError.invalidMesh("State vector length does not match element count.")
        }

        let elementCount = preparedMesh.elements.count
        var forces = Array(repeating: SIMD3<Float>.zero, count: elementCount * 4)
        var states = Array(repeating: ElementState.zero, count: elementCount)

        for elementIndex in 0..<elementCount {
            let geometry = preparedMesh.elements[elementIndex]
            let reference = gatherElementNodalVectors(geometry: geometry, values: preparedMesh.nodes)
            let localDisplacements = gatherElementNodalVectors(geometry: geometry, values: displacements)

            let response = computeLocalElementResponse(
                geometry: geometry,
                reference: reference,
                displacements: localDisplacements,
                previousState: previousStates[elementIndex],
                material: material
            )

            let base = elementIndex * 4
            forces[base] = response.forces.n0
            forces[base + 1] = response.forces.n1
            forces[base + 2] = response.forces.n2
            forces[base + 3] = response.forces.n3
            states[elementIndex] = response.updatedState
        }

        return ElementEvaluation(elementNodeForces: forces, trialStates: states)
    }
}
