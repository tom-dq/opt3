import Foundation
import simd

struct ElementNodalVectors {
    var n0: SIMD3<Float>
    var n1: SIMD3<Float>
    var n2: SIMD3<Float>
    var n3: SIMD3<Float>

    subscript(localIndex: Int) -> SIMD3<Float> {
        get {
            switch localIndex {
            case 0: return n0
            case 1: return n1
            case 2: return n2
            case 3: return n3
            default:
                preconditionFailure("Local node index must be in 0...3")
            }
        }
        set {
            switch localIndex {
            case 0: n0 = newValue
            case 1: n1 = newValue
            case 2: n2 = newValue
            case 3: n3 = newValue
            default:
                preconditionFailure("Local node index must be in 0...3")
            }
        }
    }
}

struct ElementLocalResponse {
    var forces: ElementNodalVectors
    var updatedState: ElementState
}

@inline(__always)
func gatherElementNodalVectors(geometry: ElementGeometry, values: [SIMD3<Float>]) -> ElementNodalVectors {
    let ids = geometry.nodeIDs
    return ElementNodalVectors(
        n0: values[Int(ids[0])],
        n1: values[Int(ids[1])],
        n2: values[Int(ids[2])],
        n3: values[Int(ids[3])]
    )
}

@inline(__always)
func flattenElementNodalVectors(_ vectors: ElementNodalVectors) -> [Float] {
    [
        vectors.n0.x, vectors.n0.y, vectors.n0.z,
        vectors.n1.x, vectors.n1.y, vectors.n1.z,
        vectors.n2.x, vectors.n2.y, vectors.n2.z,
        vectors.n3.x, vectors.n3.y, vectors.n3.z,
    ]
}

@inline(__always)
func computeLocalElementResponse(
    geometry: ElementGeometry,
    reference: ElementNodalVectors,
    displacements: ElementNodalVectors,
    previousState: ElementState,
    material: MaterialParameters,
    regularization: Float = 1e-4
) -> ElementLocalResponse {
    let x0 = reference.n0 + displacements.n0
    let x1 = reference.n1 + displacements.n1
    let x2 = reference.n2 + displacements.n2
    let x3 = reference.n3 + displacements.n3

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

    return ElementLocalResponse(
        forces: ElementNodalVectors(
            n0: geometry.volume * (firstPiola * g0),
            n1: geometry.volume * (firstPiola * g1),
            n2: geometry.volume * (firstPiola * g2),
            n3: geometry.volume * (firstPiola * g3)
        ),
        updatedState: ElementState(equivalentPlasticStrain: equivalentPlasticStrain, damage: damage)
    )
}
