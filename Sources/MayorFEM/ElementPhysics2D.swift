import Foundation
import simd

struct ElementNodalVectors2D {
    var n0: SIMD2<Float>
    var n1: SIMD2<Float>
    var n2: SIMD2<Float>
    var n3: SIMD2<Float>

    subscript(localIndex: Int) -> SIMD2<Float> {
        get {
            switch localIndex {
            case 0: return n0
            case 1: return n1
            case 2: return n2
            case 3: return n3
            default:
                preconditionFailure("Local quad node index must be in 0...3")
            }
        }
        set {
            switch localIndex {
            case 0: n0 = newValue
            case 1: n1 = newValue
            case 2: n2 = newValue
            case 3: n3 = newValue
            default:
                preconditionFailure("Local quad node index must be in 0...3")
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

struct ElementLinearization2D {
    var nodeIDs: SIMD4<UInt32>
    var localResidual: [Float]
    var localTangent: [Float]
    var updatedState: ElementState
    var vonMisesStress: Float
    var strainEnergy: Float
}

@inline(__always)
private func outerProduct2(_ a: SIMD2<Float>, _ b: SIMD2<Float>) -> simd_float2x2 {
    simd_float2x2(columns: (a * b.x, a * b.y))
}

@inline(__always)
func gatherElementNodalVectors2D(geometry: Quad4Geometry2D, values: [SIMD2<Float>]) -> ElementNodalVectors2D {
    let ids = geometry.nodeIDs
    return ElementNodalVectors2D(
        n0: values[Int(ids[0])],
        n1: values[Int(ids[1])],
        n2: values[Int(ids[2])],
        n3: values[Int(ids[3])]
    )
}

@inline(__always)
func flattenElementNodalVectors2D(_ vectors: ElementNodalVectors2D) -> [Float] {
    [
        vectors.n0.x, vectors.n0.y,
        vectors.n1.x, vectors.n1.y,
        vectors.n2.x, vectors.n2.y,
        vectors.n3.x, vectors.n3.y,
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

private struct QuadIntegrationPoint {
    var xi: Float
    var eta: Float
    var weight: Float
}

private func integrationPoints2D(scheme: IntegrationScheme2D) -> [QuadIntegrationPoint] {
    switch scheme {
    case .reduced:
        return [QuadIntegrationPoint(xi: 0, eta: 0, weight: 4)]
    case .full:
        let a = Float(1.0 / sqrt(3.0))
        return [
            QuadIntegrationPoint(xi: -a, eta: -a, weight: 1),
            QuadIntegrationPoint(xi: a, eta: -a, weight: 1),
            QuadIntegrationPoint(xi: a, eta: a, weight: 1),
            QuadIntegrationPoint(xi: -a, eta: a, weight: 1),
        ]
    }
}

@inline(__always)
private func shapeDerivativesQuad4(xi: Float, eta: Float) -> (dXi: [Float], dEta: [Float]) {
    let dXi: [Float] = [
        -0.25 * (1 - eta),
         0.25 * (1 - eta),
         0.25 * (1 + eta),
        -0.25 * (1 + eta),
    ]

    let dEta: [Float] = [
        -0.25 * (1 - xi),
        -0.25 * (1 + xi),
         0.25 * (1 + xi),
         0.25 * (1 - xi),
    ]

    return (dXi: dXi, dEta: dEta)
}

@inline(__always)
private func materialPointUpdate2D(
    deformationGradient2D: simd_float2x2,
    previousState: ElementState,
    material: MaterialParameters,
    regularization: Float
) -> (firstPiola2D: simd_float2x2, updatedState: ElementState, vonMises: Float) {
    var F2 = deformationGradient2D
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
    let firstPiola2D = simd_float2x2(columns: (
        SIMD2<Float>(firstPiola3D.columns.0.x, firstPiola3D.columns.0.y),
        SIMD2<Float>(firstPiola3D.columns.1.x, firstPiola3D.columns.1.y)
    ))

    return (
        firstPiola2D: firstPiola2D,
        updatedState: ElementState(equivalentPlasticStrain: equivalentPlasticStrain, damage: damage),
        vonMises: vonMisesStress
    )
}

@inline(__always)
func computeLocalElementResponse2D(
    geometry: Quad4Geometry2D,
    reference: ElementNodalVectors2D,
    displacements: ElementNodalVectors2D,
    previousState: ElementState,
    material: MaterialParameters,
    thickness: Float,
    integrationScheme: IntegrationScheme2D,
    regularization: Float = 1e-5
) -> ElementLocalResponse2D {
    let referenceNodes = [reference.n0, reference.n1, reference.n2, reference.n3]
    let displacementNodes = [displacements.n0, displacements.n1, displacements.n2, displacements.n3]
    let currentNodes = [
        reference.n0 + displacements.n0,
        reference.n1 + displacements.n1,
        reference.n2 + displacements.n2,
        reference.n3 + displacements.n3,
    ]

    let points = integrationPoints2D(scheme: integrationScheme)
    var nodalForces = Array(repeating: SIMD2<Float>.zero, count: 4)

    var weightedEqp: Float = 0
    var weightedDamage: Float = 0
    var weightedVonMises: Float = 0
    var totalWeight: Float = 0

    for point in points {
        let derivatives = shapeDerivativesQuad4(xi: point.xi, eta: point.eta)

        var dX_dXi = SIMD2<Float>.zero
        var dX_dEta = SIMD2<Float>.zero
        for i in 0..<4 {
            dX_dXi += derivatives.dXi[i] * referenceNodes[i]
            dX_dEta += derivatives.dEta[i] * referenceNodes[i]
        }

        var jacobian = simd_float2x2(columns: (dX_dXi, dX_dEta))
        var detJ = simd_determinant(jacobian)
        if detJ < regularization {
            jacobian += (regularization - detJ) * matrix_identity_float2x2
            detJ = simd_determinant(jacobian)
        }
        detJ = max(detJ, regularization)

        let inverseTransposeJ = simd_transpose(simd_inverse(jacobian))

        var gradients: [SIMD2<Float>] = []
        gradients.reserveCapacity(4)
        for i in 0..<4 {
            gradients.append(inverseTransposeJ * SIMD2<Float>(derivatives.dXi[i], derivatives.dEta[i]))
        }

        var deformationGradient2D = simd_float2x2(columns: (SIMD2<Float>.zero, SIMD2<Float>.zero))
        for i in 0..<4 {
            deformationGradient2D += outerProduct2(currentNodes[i], gradients[i])
        }

        let update = materialPointUpdate2D(
            deformationGradient2D: deformationGradient2D,
            previousState: previousState,
            material: material,
            regularization: regularization
        )

        let weight = thickness * detJ * point.weight
        for i in 0..<4 {
            nodalForces[i] += weight * (update.firstPiola2D * gradients[i])
        }

        weightedEqp += weight * update.updatedState.equivalentPlasticStrain
        weightedDamage += weight * update.updatedState.damage
        weightedVonMises += weight * update.vonMises
        totalWeight += weight
    }

    let invWeight = totalWeight > 0 ? 1.0 / totalWeight : 0
    let updatedState = ElementState(
        equivalentPlasticStrain: weightedEqp * invWeight,
        damage: weightedDamage * invWeight
    )
    let vonMises = weightedVonMises * invWeight

    var strainEnergy: Float = 0
    for i in 0..<4 {
        strainEnergy += simd_dot(nodalForces[i], displacementNodes[i])
    }
    strainEnergy = abs(0.5 * strainEnergy)

    return ElementLocalResponse2D(
        forces: ElementNodalVectors2D(
            n0: nodalForces[0],
            n1: nodalForces[1],
            n2: nodalForces[2],
            n3: nodalForces[3]
        ),
        updatedState: updatedState,
        vonMisesStress: vonMises,
        strainEnergy: strainEnergy
    )
}

final class CPUElementLinearizer2D {
    private let preparedMesh: PreparedMesh2D
    private let material: MaterialParameters
    private let thickness: Float
    private let integrationScheme: IntegrationScheme2D

    init(
        preparedMesh: PreparedMesh2D,
        material: MaterialParameters,
        thickness: Float,
        integrationScheme: IntegrationScheme2D
    ) {
        self.preparedMesh = preparedMesh
        self.material = material
        self.thickness = thickness
        self.integrationScheme = integrationScheme
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
            thickness: thickness,
            integrationScheme: integrationScheme
        )
        let baseResidual = flattenElementNodalVectors2D(baseResponse.forces)
        var localTangent = Array(repeating: Float(0), count: 8 * 8)

        for column in 0..<8 {
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
                thickness: thickness,
                integrationScheme: integrationScheme
            )
            let perturbedResidual = flattenElementNodalVectors2D(perturbedResponse.forces)

            for row in 0..<8 {
                localTangent[row * 8 + column] =
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
