import Foundation
import simd

struct ElementNodalVectors2D {
    var n0: SIMD2<Float>
    var n1: SIMD2<Float>
    var n2: SIMD2<Float>

    subscript(localIndex: Int) -> SIMD2<Float> {
        get {
            switch localIndex {
            case 0: return n0
            case 1: return n1
            case 2: return n2
            default:
                preconditionFailure("Local node index must be in 0...2")
            }
        }
        set {
            switch localIndex {
            case 0: n0 = newValue
            case 1: n1 = newValue
            case 2: n2 = newValue
            default:
                preconditionFailure("Local node index must be in 0...2")
            }
        }
    }
}

struct ElementLocalResponse2D {
    var forces: ElementNodalVectors2D
    var updatedState: ElementState
    var vonMisesStress: Float
    var strainEnergy: Float
}

@inline(__always)
private func outerProduct2(_ a: SIMD2<Float>, _ b: SIMD2<Float>) -> simd_float2x2 {
    simd_float2x2(columns: (a * b.x, a * b.y))
}

@inline(__always)
func gatherElementNodalVectors2D(geometry: Tri3Geometry2D, values: [SIMD2<Float>]) -> ElementNodalVectors2D {
    let ids = geometry.nodeIDs
    return ElementNodalVectors2D(
        n0: values[Int(ids[0])],
        n1: values[Int(ids[1])],
        n2: values[Int(ids[2])]
    )
}

@inline(__always)
func flattenElementNodalVectors2D(_ vectors: ElementNodalVectors2D) -> [Float] {
    [
        vectors.n0.x, vectors.n0.y,
        vectors.n1.x, vectors.n1.y,
        vectors.n2.x, vectors.n2.y,
    ]
}

@inline(__always)
func dofVectorToNodalDisplacements2D(_ dofs: [Float], nodeCount: Int) -> [SIMD2<Float>] {
    var displacements = Array(repeating: SIMD2<Float>.zero, count: nodeCount)
    for node in 0..<nodeCount {
        let base = node * 2
        displacements[node] = SIMD2<Float>(dofs[base], dofs[base + 1])
    }
    return displacements
}

@inline(__always)
private func embedPlaneStrain(_ F2: simd_float2x2) -> simd_float3x3 {
    simd_float3x3(columns: (
        SIMD3<Float>(F2.columns.0.x, F2.columns.0.y, 0),
        SIMD3<Float>(F2.columns.1.x, F2.columns.1.y, 0),
        SIMD3<Float>(0, 0, 1)
    ))
}

@inline(__always)
func computeLocalElementResponse2D(
    geometry: Tri3Geometry2D,
    reference: ElementNodalVectors2D,
    displacements: ElementNodalVectors2D,
    previousState: ElementState,
    material: MaterialParameters,
    thickness: Float,
    regularization: Float = 1e-5
) -> ElementLocalResponse2D {
    let x0 = reference.n0 + displacements.n0
    let x1 = reference.n1 + displacements.n1
    let x2 = reference.n2 + displacements.n2

    let g0 = geometry.gradN0
    let g1 = geometry.gradN1
    let g2 = geometry.gradN2

    var F2 = outerProduct2(x0, g0)
    F2 += outerProduct2(x1, g1)
    F2 += outerProduct2(x2, g2)

    var detF2 = simd_determinant(F2)
    if detF2 < regularization {
        F2 += (regularization - detF2) * matrix_identity_float2x2
        detF2 = simd_determinant(F2)
    }
    detF2 = max(detF2, regularization)

    let F = embedPlaneStrain(F2)

    let leftCauchyGreen = F * simd_transpose(F)
    let identity = matrix_identity_float3x3
    let trialKirchhoff = material.mu * (leftCauchyGreen - identity) + material.lambda * log(detF2) * identity

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
    let deviatoricKirchhoff = softening * deviatoric
    let vonMisesStress = sqrt(max(0, 1.5 * doubleContraction(deviatoricKirchhoff)))

    let inverseTransposeF = simd_transpose(simd_inverse(F))
    let firstPiola3D = kirchhoff * inverseTransposeF
    let firstPiola = simd_float2x2(columns: (
        SIMD2<Float>(firstPiola3D.columns.0.x, firstPiola3D.columns.0.y),
        SIMD2<Float>(firstPiola3D.columns.1.x, firstPiola3D.columns.1.y)
    ))

    let scale = thickness * geometry.area

    return ElementLocalResponse2D(
        forces: ElementNodalVectors2D(
            n0: scale * (firstPiola * g0),
            n1: scale * (firstPiola * g1),
            n2: scale * (firstPiola * g2)
        ),
        updatedState: ElementState(equivalentPlasticStrain: equivalentPlasticStrain, damage: damage),
        vonMisesStress: vonMisesStress,
        strainEnergy: abs(
            0.5 * (
                simd_dot(scale * (firstPiola * g0), displacements.n0) +
                simd_dot(scale * (firstPiola * g1), displacements.n1) +
                simd_dot(scale * (firstPiola * g2), displacements.n2)
            )
        )
    )
}

struct ElementLinearization2D {
    var nodeIDs: SIMD3<UInt32>
    var localResidual: [Float]
    var localTangent: [Float]
    var updatedState: ElementState
    var vonMisesStress: Float
    var strainEnergy: Float
}

final class CPUElementLinearizer2D {
    private let preparedMesh: PreparedMesh2D
    private let material: MaterialParameters
    private let thickness: Float

    init(preparedMesh: PreparedMesh2D, material: MaterialParameters, thickness: Float) {
        self.preparedMesh = preparedMesh
        self.material = material
        self.thickness = thickness
    }

    func linearize(
        elementIndex: Int,
        displacements: [SIMD2<Float>],
        previousState: ElementState,
        finiteDifferenceStep: Float
    ) -> ElementLinearization2D {
        let geometry = preparedMesh.elements[elementIndex]
        let reference = gatherElementNodalVectors2D(geometry: geometry, values: preparedMesh.nodes)
        let localDisplacements = gatherElementNodalVectors2D(geometry: geometry, values: displacements)

        let baseResponse = computeLocalElementResponse2D(
            geometry: geometry,
            reference: reference,
            displacements: localDisplacements,
            previousState: previousState,
            material: material,
            thickness: thickness
        )
        let baseResidual = flattenElementNodalVectors2D(baseResponse.forces)
        var localTangent = Array(repeating: Float(0), count: 6 * 6)

        for column in 0..<6 {
            let node = column / 2
            let component = column % 2

            var perturbedDisplacements = localDisplacements
            var vector = perturbedDisplacements[node]
            let perturbation = finiteDifferenceStep * max(1.0, abs(vector[component]))
            vector[component] += perturbation
            perturbedDisplacements[node] = vector

            let perturbedResponse = computeLocalElementResponse2D(
                geometry: geometry,
                reference: reference,
                displacements: perturbedDisplacements,
                previousState: previousState,
                material: material,
                thickness: thickness
            )
            let perturbedResidual = flattenElementNodalVectors2D(perturbedResponse.forces)

            for row in 0..<6 {
                localTangent[row * 6 + column] =
                    (perturbedResidual[row] - baseResidual[row]) / perturbation
            }
        }

        return ElementLinearization2D(
            nodeIDs: geometry.nodeIDs,
            localResidual: baseResidual,
            localTangent: localTangent,
            updatedState: baseResponse.updatedState,
            vonMisesStress: baseResponse.vonMisesStress,
            strainEnergy: baseResponse.strainEnergy
        )
    }
}
