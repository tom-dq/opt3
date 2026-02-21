import Foundation
import simd

struct ElementEvaluation2D {
    var elementNodeForces: [SIMD2<Float>]
    var trialStates: [ElementState]
    var elementVonMises: [Float]
    var elementEnergies: [Float]
}

protocol ElementEvaluator2D {
    func evaluate(
        displacements: [SIMD2<Float>],
        previousStates: [ElementState],
        densities: [Float],
        densityPenalty: Float,
        minimumDensity: Float
    ) throws -> ElementEvaluation2D
}

@inline(__always)
func densityScaleFactor(_ density: Float, penalty: Float, minimumDensity: Float) -> Float {
    let clamped = max(minimumDensity, min(1.0, density))
    return pow(clamped, max(1.0, penalty))
}

final class CPUElementEvaluator2D: ElementEvaluator2D {
    private let preparedMesh: PreparedMesh2D
    private let material: MaterialParameters
    private let thickness: Float

    init(preparedMesh: PreparedMesh2D, material: MaterialParameters, thickness: Float) {
        self.preparedMesh = preparedMesh
        self.material = material
        self.thickness = thickness
    }

    func evaluate(
        displacements: [SIMD2<Float>],
        previousStates: [ElementState],
        densities: [Float],
        densityPenalty: Float,
        minimumDensity: Float
    ) throws -> ElementEvaluation2D {
        if displacements.count != preparedMesh.nodes.count {
            throw FEMError.invalidMesh("2D displacement vector length does not match node count.")
        }
        if previousStates.count != preparedMesh.elements.count {
            throw FEMError.invalidMesh("2D state vector length does not match element count.")
        }
        if densities.count != preparedMesh.elements.count {
            throw FEMError.invalidMesh("2D density vector length does not match element count.")
        }

        let elementCount = preparedMesh.elements.count
        var forces = Array(repeating: SIMD2<Float>.zero, count: elementCount * 3)
        var states = Array(repeating: ElementState.zero, count: elementCount)
        var vonMises = Array(repeating: Float(0), count: elementCount)
        var energies = Array(repeating: Float(0), count: elementCount)

        for elementIndex in 0..<elementCount {
            let geometry = preparedMesh.elements[elementIndex]
            let reference = gatherElementNodalVectors2D(geometry: geometry, values: preparedMesh.nodes)
            let localDisplacements = gatherElementNodalVectors2D(geometry: geometry, values: displacements)
            let previousState = previousStates[elementIndex]
            let response = computeLocalElementResponse2D(
                geometry: geometry,
                reference: reference,
                displacements: localDisplacements,
                previousState: previousState,
                material: material,
                thickness: thickness
            )

            let scale = densityScaleFactor(densities[elementIndex], penalty: densityPenalty, minimumDensity: minimumDensity)

            let base = elementIndex * 3
            forces[base] = scale * response.forces.n0
            forces[base + 1] = scale * response.forces.n1
            forces[base + 2] = scale * response.forces.n2

            let blendedEqp = previousState.equivalentPlasticStrain +
                scale * (response.updatedState.equivalentPlasticStrain - previousState.equivalentPlasticStrain)
            let blendedDamage = previousState.damage +
                scale * (response.updatedState.damage - previousState.damage)
            states[elementIndex] = ElementState(equivalentPlasticStrain: blendedEqp, damage: blendedDamage)
            vonMises[elementIndex] = scale * response.vonMisesStress
            energies[elementIndex] = scale * response.strainEnergy
        }

        return ElementEvaluation2D(
            elementNodeForces: forces,
            trialStates: states,
            elementVonMises: vonMises,
            elementEnergies: energies
        )
    }
}
