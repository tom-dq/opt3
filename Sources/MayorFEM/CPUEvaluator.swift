import Foundation
import simd

final class CPUEvaluator: ElementEvaluator {
    private let preparedMesh: PreparedMesh
    private let material: MaterialParameters
    private let regularization: Float = 1e-4

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
            let response = computeElementResponse(
                geometry: preparedMesh.elements[elementIndex],
                displacements: displacements,
                previousState: previousStates[elementIndex]
            )

            let base = elementIndex * 4
            forces[base] = response.forces[0]
            forces[base + 1] = response.forces[1]
            forces[base + 2] = response.forces[2]
            forces[base + 3] = response.forces[3]
            states[elementIndex] = response.updatedState
        }

        return ElementEvaluation(elementNodeForces: forces, trialStates: states)
    }

    private func computeElementResponse(
        geometry: ElementGeometry,
        displacements: [SIMD3<Float>],
        previousState: ElementState
    ) -> (forces: [SIMD3<Float>], updatedState: ElementState) {
        let ids = geometry.nodeIDs

        let x0 = preparedMesh.nodes[Int(ids[0])] + displacements[Int(ids[0])]
        let x1 = preparedMesh.nodes[Int(ids[1])] + displacements[Int(ids[1])]
        let x2 = preparedMesh.nodes[Int(ids[2])] + displacements[Int(ids[2])]
        let x3 = preparedMesh.nodes[Int(ids[3])] + displacements[Int(ids[3])]

        let g0 = geometry.gradient(localIndex: 0)
        let g1 = geometry.gradient(localIndex: 1)
        let g2 = geometry.gradient(localIndex: 2)
        let g3 = geometry.gradient(localIndex: 3)

        var deformationGradient = outerProduct(x0, g0)
        deformationGradient += outerProduct(x1, g1)
        deformationGradient += outerProduct(x2, g2)
        deformationGradient += outerProduct(x3, g3)

        let identity = matrix_identity_float3x3

        var determinantF = simd_determinant(deformationGradient)
        if determinantF < regularization {
            deformationGradient += (regularization - determinantF) * identity
            determinantF = simd_determinant(deformationGradient)
        }
        determinantF = max(determinantF, regularization)

        let leftCauchyGreen = deformationGradient * simd_transpose(deformationGradient)
        let trialKirchhoff = material.mu * (leftCauchyGreen - identity) + material.lambda * log(determinantF) * identity

        let meanStress = trace(trialKirchhoff) / 3.0
        var deviatoric = trialKirchhoff - meanStress * identity

        let qTrial = sqrt(max(1e-12, 1.5 * doubleContraction(deviatoric)))
        let flowStress = material.yieldStress + material.hardeningModulus * previousState.equivalentPlasticStrain

        var equivalentPlasticStrain = previousState.equivalentPlasticStrain
        if qTrial > flowStress {
            let denominator = 3.0 * material.mu + material.hardeningModulus + 1e-6
            let deltaGamma = (qTrial - flowStress) / denominator
            let scale = max(0, 1.0 - (3.0 * material.mu * deltaGamma) / (qTrial + 1e-8))
            deviatoric *= scale
            equivalentPlasticStrain += sqrt(2.0 / 3.0) * deltaGamma
        }

        var damage = previousState.damage
        if equivalentPlasticStrain > material.damageOnset {
            let threshold = max(previousState.equivalentPlasticStrain, material.damageOnset)
            let deltaEquivalentPlastic = max(0, equivalentPlasticStrain - threshold)
            let deltaDamage = deltaEquivalentPlastic / max(material.damageSlope, 1e-6)
            damage = min(material.damageCap, previousState.damage + deltaDamage)
        }

        let softening = max(0, 1.0 - damage)
        let kirchhoff = softening * (deviatoric + meanStress * identity)

        let inverseTransposeF = simd_transpose(simd_inverse(deformationGradient))
        let firstPiola = kirchhoff * inverseTransposeF

        let force0 = geometry.volume * (firstPiola * g0)
        let force1 = geometry.volume * (firstPiola * g1)
        let force2 = geometry.volume * (firstPiola * g2)
        let force3 = geometry.volume * (firstPiola * g3)

        return (
            forces: [force0, force1, force2, force3],
            updatedState: ElementState(equivalentPlasticStrain: equivalentPlasticStrain, damage: damage)
        )
    }
}
