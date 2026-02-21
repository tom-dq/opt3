import Foundation
import simd

private struct GlobalEvaluation {
    var residual: [Float]
    var trialStates: [ElementState]
}

public final class NonlinearFEMSolver {
    private let problem: FEMProblem
    private let preparedMesh: PreparedMesh
    private let evaluator: ElementEvaluator
    private let backendName: String
    private let dofCount: Int
    private let freeDOFs: [Int]
    private let prescribedFinalValues: [Int: Float]

    public init(problem: FEMProblem, backendChoice: ComputeBackendChoice = .auto) throws {
        self.problem = problem
        self.preparedMesh = try problem.mesh.prepare()
        self.dofCount = preparedMesh.nodes.count * 3

        let prescribedFinalValues = try NonlinearFEMSolver.buildPrescribedMap(
            prescribedDisplacements: problem.prescribedDisplacements,
            dofCount: dofCount
        )
        self.prescribedFinalValues = prescribedFinalValues

        let prescribedSet = Set(prescribedFinalValues.keys)
        self.freeDOFs = (0..<dofCount).filter { !prescribedSet.contains($0) }

        if freeDOFs.isEmpty {
            throw FEMError.noFreeDOFs
        }

        switch backendChoice {
        case .cpu:
            self.evaluator = CPUEvaluator(preparedMesh: preparedMesh, material: problem.material)
            self.backendName = "cpu"
        case .metal:
            guard let metalEvaluator = try NonlinearFEMSolver.makeMetalEvaluator(preparedMesh: preparedMesh, material: problem.material) else {
                throw FEMError.backendUnavailable("Metal backend requested but unavailable.")
            }
            self.evaluator = metalEvaluator
            self.backendName = "metal"
        case .auto:
            if let metalEvaluator = try NonlinearFEMSolver.makeMetalEvaluator(preparedMesh: preparedMesh, material: problem.material) {
                self.evaluator = metalEvaluator
                self.backendName = "metal"
            } else {
                self.evaluator = CPUEvaluator(preparedMesh: preparedMesh, material: problem.material)
                self.backendName = "cpu"
            }
        }
    }

    public func solve() throws -> SolveResult {
        var dofs = Array(repeating: Float(0), count: dofCount)
        var states = Array(repeating: ElementState.zero, count: preparedMesh.elements.count)
        var stepHistory: [StepResult] = []
        var finalReactions = Array(repeating: Float(0), count: dofCount)

        for step in 1...problem.controls.loadSteps {
            let loadFactor = Float(step) / Float(problem.controls.loadSteps)
            let prescribedAtStep = prescribedValues(for: loadFactor)
            applyPrescribedValues(&dofs, prescribed: prescribedAtStep)

            let baseStates = states
            var convergedEvaluation: GlobalEvaluation?
            var iterations = 0
            var residualNorm: Float = .greatestFiniteMagnitude

            for iteration in 1...problem.controls.maxNewtonIterations {
                iterations = iteration
                let evaluation = try evaluateGlobal(dofs: dofs, previousStates: baseStates)
                residualNorm = l2Norm(evaluation.residual, indices: freeDOFs)

                if residualNorm <= problem.controls.residualTolerance {
                    convergedEvaluation = evaluation
                    break
                }

                let jacobian = try buildNumericalJacobian(
                    dofs: dofs,
                    baseResidual: evaluation.residual,
                    prescribedAtStep: prescribedAtStep,
                    previousStates: baseStates
                )
                let rhs = freeDOFs.map { -evaluation.residual[$0] }
                let increment = try solveDenseLinearSystem(coefficients: jacobian, rhs: rhs)

                var accepted = false
                var alpha: Float = 1.0
                while alpha >= problem.controls.lineSearchFloor {
                    var candidate = dofs
                    for (index, dofIndex) in freeDOFs.enumerated() {
                        candidate[dofIndex] += alpha * increment[index]
                    }
                    applyPrescribedValues(&candidate, prescribed: prescribedAtStep)

                    let candidateEvaluation = try evaluateGlobal(dofs: candidate, previousStates: baseStates)
                    let candidateNorm = l2Norm(candidateEvaluation.residual, indices: freeDOFs)
                    if candidateNorm < residualNorm {
                        dofs = candidate
                        accepted = true
                        break
                    }

                    alpha *= 0.5
                }

                if !accepted {
                    for (index, dofIndex) in freeDOFs.enumerated() {
                        dofs[dofIndex] += 0.1 * increment[index]
                    }
                    applyPrescribedValues(&dofs, prescribed: prescribedAtStep)
                }
            }

            if convergedEvaluation == nil {
                let terminalEvaluation = try evaluateGlobal(dofs: dofs, previousStates: baseStates)
                residualNorm = l2Norm(terminalEvaluation.residual, indices: freeDOFs)
                throw FEMError.newtonFailed(step: step, residualNorm: residualNorm)
            }

            states = convergedEvaluation!.trialStates
            finalReactions = convergedEvaluation!.residual

            let maxEquivalentPlastic = states.map(\.equivalentPlasticStrain).max() ?? 0
            let maxDamage = states.map(\.damage).max() ?? 0
            stepHistory.append(
                StepResult(
                    step: step,
                    loadFactor: loadFactor,
                    iterations: iterations,
                    residualNorm: residualNorm,
                    maxEquivalentPlasticStrain: maxEquivalentPlastic,
                    maxDamage: maxDamage
                )
            )
        }

        return SolveResult(
            displacements: dofVectorToNodalDisplacements(dofs, nodeCount: preparedMesh.nodes.count),
            reactions: finalReactions,
            elementStates: states,
            stepHistory: stepHistory,
            backendName: backendName,
            converged: true
        )
    }

    private func evaluateGlobal(dofs: [Float], previousStates: [ElementState]) throws -> GlobalEvaluation {
        let nodalDisplacements = dofVectorToNodalDisplacements(dofs, nodeCount: preparedMesh.nodes.count)
        let elementEvaluation = try evaluator.evaluate(
            displacements: nodalDisplacements,
            previousStates: previousStates
        )

        var residual = Array(repeating: Float(0), count: dofCount)

        for (elementIndex, element) in preparedMesh.elements.enumerated() {
            let forceBase = elementIndex * 4
            for localIndex in 0..<4 {
                let node = Int(element.nodeIDs[localIndex])
                let force = elementEvaluation.elementNodeForces[forceBase + localIndex]
                let dofBase = node * 3
                residual[dofBase] += force.x
                residual[dofBase + 1] += force.y
                residual[dofBase + 2] += force.z
            }
        }

        return GlobalEvaluation(residual: residual, trialStates: elementEvaluation.trialStates)
    }

    private func buildNumericalJacobian(
        dofs: [Float],
        baseResidual: [Float],
        prescribedAtStep: [Int: Float],
        previousStates: [ElementState]
    ) throws -> [Float] {
        let dimension = freeDOFs.count
        var jacobian = Array(repeating: Float(0), count: dimension * dimension)

        for column in 0..<dimension {
            let dofIndex = freeDOFs[column]
            var perturbed = dofs
            let perturbation = problem.controls.finiteDifferenceStep * max(1, abs(dofs[dofIndex]))
            perturbed[dofIndex] += perturbation
            applyPrescribedValues(&perturbed, prescribed: prescribedAtStep)

            let perturbedResidual = try evaluateGlobal(dofs: perturbed, previousStates: previousStates).residual

            for row in 0..<dimension {
                let residualIndex = freeDOFs[row]
                jacobian[row * dimension + column] =
                    (perturbedResidual[residualIndex] - baseResidual[residualIndex]) / perturbation
            }
        }

        return jacobian
    }

    private func prescribedValues(for loadFactor: Float) -> [Int: Float] {
        var values: [Int: Float] = [:]
        values.reserveCapacity(prescribedFinalValues.count)
        for (dof, finalValue) in prescribedFinalValues {
            values[dof] = finalValue * loadFactor
        }
        return values
    }

    private static func buildPrescribedMap(
        prescribedDisplacements: [PrescribedDisplacement],
        dofCount: Int
    ) throws -> [Int: Float] {
        var map: [Int: Float] = [:]
        for constraint in prescribedDisplacements {
            if !(0..<3).contains(constraint.component) {
                throw FEMError.invalidBoundaryCondition(
                    "Node \(constraint.node) uses invalid displacement component \(constraint.component)."
                )
            }

            let dof = constraint.dof
            if dof < 0 || dof >= dofCount {
                throw FEMError.invalidBoundaryCondition("DOF \(dof) is out of range.")
            }

            if let existing = map[dof], abs(existing - constraint.value) > 1e-8 {
                throw FEMError.invalidBoundaryCondition("Conflicting prescribed displacement values for DOF \(dof).")
            }

            map[dof] = constraint.value
        }
        return map
    }

    private static func makeMetalEvaluator(
        preparedMesh: PreparedMesh,
        material: MaterialParameters
    ) throws -> ElementEvaluator? {
        #if canImport(Metal)
        return try MetalElementEvaluator(preparedMesh: preparedMesh, material: material)
        #else
        _ = preparedMesh
        _ = material
        return nil
        #endif
    }
}
