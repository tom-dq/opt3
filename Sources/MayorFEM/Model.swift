import Foundation
import simd

public struct Mesh {
    public var nodes: [SIMD3<Float>]
    public var elements: [SIMD4<UInt32>]

    public init(nodes: [SIMD3<Float>], elements: [SIMD4<UInt32>]) {
        self.nodes = nodes
        self.elements = elements
    }

    public static func fiveTetraBlock(length: Float = 1.0, width: Float = 0.2, height: Float = 0.2) -> Mesh {
        let nodes: [SIMD3<Float>] = [
            SIMD3<Float>(0, 0, 0),
            SIMD3<Float>(length, 0, 0),
            SIMD3<Float>(length, width, 0),
            SIMD3<Float>(0, width, 0),
            SIMD3<Float>(0, 0, height),
            SIMD3<Float>(length, 0, height),
            SIMD3<Float>(length, width, height),
            SIMD3<Float>(0, width, height),
        ]

        let elements: [SIMD4<UInt32>] = [
            SIMD4<UInt32>(0, 1, 3, 4),
            SIMD4<UInt32>(1, 2, 3, 6),
            SIMD4<UInt32>(1, 4, 5, 6),
            SIMD4<UInt32>(3, 4, 6, 7),
            SIMD4<UInt32>(1, 3, 4, 6),
        ]

        return Mesh(nodes: nodes, elements: elements)
    }
}

public struct MaterialParameters {
    public var youngsModulus: Float
    public var poissonRatio: Float
    public var yieldStress: Float
    public var hardeningModulus: Float
    public var damageOnset: Float
    public var damageSlope: Float
    public var damageCap: Float

    public init(
        youngsModulus: Float,
        poissonRatio: Float,
        yieldStress: Float,
        hardeningModulus: Float,
        damageOnset: Float,
        damageSlope: Float,
        damageCap: Float
    ) {
        self.youngsModulus = youngsModulus
        self.poissonRatio = poissonRatio
        self.yieldStress = yieldStress
        self.hardeningModulus = hardeningModulus
        self.damageOnset = damageOnset
        self.damageSlope = damageSlope
        self.damageCap = damageCap
    }

    var lambda: Float {
        let denominator = (1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)
        return (youngsModulus * poissonRatio) / max(denominator, 1e-6)
    }

    var mu: Float {
        youngsModulus / (2.0 * (1.0 + poissonRatio))
    }
}

public struct PrescribedDisplacement {
    public var node: Int
    public var component: Int
    public var value: Float

    public init(node: Int, component: Int, value: Float) {
        self.node = node
        self.component = component
        self.value = value
    }

    public var dof: Int {
        node * 3 + component
    }
}

public struct SolverControls {
    public var loadSteps: Int
    public var maxNewtonIterations: Int
    public var residualTolerance: Float
    public var finiteDifferenceStep: Float
    public var lineSearchFloor: Float

    public init(
        loadSteps: Int = 12,
        maxNewtonIterations: Int = 20,
        residualTolerance: Float = 1e-4,
        finiteDifferenceStep: Float = 1e-4,
        lineSearchFloor: Float = 1.0 / 128.0
    ) {
        self.loadSteps = loadSteps
        self.maxNewtonIterations = maxNewtonIterations
        self.residualTolerance = residualTolerance
        self.finiteDifferenceStep = finiteDifferenceStep
        self.lineSearchFloor = lineSearchFloor
    }
}

public struct FEMProblem {
    public var mesh: Mesh
    public var material: MaterialParameters
    public var prescribedDisplacements: [PrescribedDisplacement]
    public var controls: SolverControls

    public init(
        mesh: Mesh,
        material: MaterialParameters,
        prescribedDisplacements: [PrescribedDisplacement],
        controls: SolverControls = SolverControls()
    ) {
        self.mesh = mesh
        self.material = material
        self.prescribedDisplacements = prescribedDisplacements
        self.controls = controls
    }
}

public struct ElementState {
    public var equivalentPlasticStrain: Float
    public var damage: Float

    public init(equivalentPlasticStrain: Float = 0, damage: Float = 0) {
        self.equivalentPlasticStrain = equivalentPlasticStrain
        self.damage = damage
    }

    public static var zero: ElementState {
        ElementState(equivalentPlasticStrain: 0, damage: 0)
    }
}

public struct StepResult {
    public var step: Int
    public var loadFactor: Float
    public var iterations: Int
    public var residualNorm: Float
    public var maxEquivalentPlasticStrain: Float
    public var maxDamage: Float
}

public struct SolveResult {
    public var displacements: [SIMD3<Float>]
    public var reactions: [Float]
    public var elementStates: [ElementState]
    public var stepHistory: [StepResult]
    public var backendName: String
    public var converged: Bool
}

public enum ComputeBackendChoice {
    case auto
    case metal
    case cpu
}

public enum FEMError: Error, CustomStringConvertible {
    case invalidMesh(String)
    case invalidBoundaryCondition(String)
    case noFreeDOFs
    case backendUnavailable(String)
    case linearSolveFailure(String)
    case newtonFailed(step: Int, residualNorm: Float)

    public var description: String {
        switch self {
        case .invalidMesh(let message):
            return "Invalid mesh: \(message)"
        case .invalidBoundaryCondition(let message):
            return "Invalid boundary condition: \(message)"
        case .noFreeDOFs:
            return "No free DOFs remain after applying displacement constraints."
        case .backendUnavailable(let message):
            return "Backend unavailable: \(message)"
        case .linearSolveFailure(let message):
            return "Linear solver failure: \(message)"
        case .newtonFailed(let step, let residualNorm):
            return "Newton solver failed at load step \(step) with residual norm \(residualNorm)."
        }
    }
}
