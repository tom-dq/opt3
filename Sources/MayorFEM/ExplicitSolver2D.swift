import Foundation
import simd

private struct ExplicitGlobalEvaluation2D {
    var residual: [Float]
    var trialStates: [ElementState]
    var elementVonMises: [Float]
    var elementEnergies: [Float]
}

public final class ExplicitFEMSolver2D {
    private let problem: FEMProblem2D
    private let explicitControls: ExplicitSolverControls2D
    private let preparedMesh: PreparedMesh2D
    private let evaluator: ElementEvaluator2D
    private let backendName: String
    private let dofCount: Int
    private let freeDOFs: [Int]
    private let prescribedFinalValues: [Int: Float]

    public init(
        problem: FEMProblem2D,
        explicitControls: ExplicitSolverControls2D = ExplicitSolverControls2D(),
        backendChoice: ComputeBackendChoice = .auto
    ) throws {
        self.problem = problem
        self.explicitControls = explicitControls
        self.preparedMesh = try problem.mesh.prepare()
        self.dofCount = preparedMesh.nodes.count * 2

        self.prescribedFinalValues = try ExplicitFEMSolver2D.buildPrescribedMap(
            prescribedDisplacements: problem.prescribedDisplacements,
            dofCount: dofCount
        )

        let prescribedSet = Set(prescribedFinalValues.keys)
        self.freeDOFs = (0..<dofCount).filter { !prescribedSet.contains($0) }
        if freeDOFs.isEmpty {
            throw FEMError.noFreeDOFs
        }

        switch backendChoice {
        case .cpu:
            self.evaluator = CPUElementEvaluator2D(
                preparedMesh: preparedMesh,
                material: problem.material,
                thickness: problem.thickness
            )
            self.backendName = "cpu-2d-explicit"
        case .metal:
            guard let metalEvaluator = try ExplicitFEMSolver2D.makeMetalEvaluator(
                preparedMesh: preparedMesh,
                material: problem.material,
                thickness: problem.thickness
            ) else {
                throw FEMError.backendUnavailable("Metal backend requested but unavailable for 2D explicit solver.")
            }
            self.evaluator = metalEvaluator
            self.backendName = "metal-2d-explicit"
        case .auto:
            if let metalEvaluator = try ExplicitFEMSolver2D.makeMetalEvaluator(
                preparedMesh: preparedMesh,
                material: problem.material,
                thickness: problem.thickness
            ) {
                self.evaluator = metalEvaluator
                self.backendName = "metal-2d-explicit"
            } else {
                self.evaluator = CPUElementEvaluator2D(
                    preparedMesh: preparedMesh,
                    material: problem.material,
                    thickness: problem.thickness
                )
                self.backendName = "cpu-2d-explicit"
            }
        }
    }

    public func solve(
        densities rawDensities: [Float]? = nil,
        minimumDensity: Float = 0.01
    ) throws -> SolveResult2D {
        let elementCount = preparedMesh.elements.count
        var densities = rawDensities ?? Array(repeating: Float(1), count: elementCount)
        if densities.count != elementCount {
            throw FEMError.invalidMesh("2D density count \(densities.count) does not match element count \(elementCount).")
        }
        let clampedMinimumDensity = max(1e-4, min(0.99, minimumDensity))
        for index in densities.indices {
            densities[index] = max(clampedMinimumDensity, min(1.0, densities[index]))
        }

        let masses = buildLumpedMasses(densities: densities, minimumDensity: clampedMinimumDensity)

        var dofs = Array(repeating: Float(0), count: dofCount)
        var states = Array(repeating: ElementState.zero, count: elementCount)
        var elementVonMises = Array(repeating: Float(0), count: elementCount)
        var elementEnergies = Array(repeating: Float(0), count: elementCount)

        var stepHistory: [StepResult2D] = []
        var stepSnapshots: [StepSnapshot2D] = []
        var finalReactions = Array(repeating: Float(0), count: dofCount)

        let dampingFactor = min(0.999, max(0.0, explicitControls.damping))
        let dt = explicitControls.timeStep
        let displacementClamp = dt * explicitControls.velocityClamp
        let substeps = explicitControls.substepsPerLoadStep
        let relaxIterations = explicitControls.relaxationIterationsPerSubstep
        let totalSubsteps = problem.controls.loadSteps * substeps

        for step in 1...problem.controls.loadSteps {
            for localSubstep in 1...substeps {
                let globalSubstep = (step - 1) * substeps + localSubstep
                let loadFactor = Float(globalSubstep) / Float(totalSubsteps)
                let prescribed = prescribedValues(for: loadFactor)

                for _ in 0..<relaxIterations {
                    applyPrescribedValues(&dofs, prescribed: prescribed)
                    let preEvaluation = try evaluateGlobal(
                        dofs: dofs,
                        previousStates: states,
                        densities: densities,
                        minimumDensity: clampedMinimumDensity
                    )

                    for dof in freeDOFs {
                        let mass = max(1e-8, masses[dof])
                        var displacementIncrement = -dt * preEvaluation.residual[dof] / mass
                        displacementIncrement *= (1.0 - dampingFactor)
                        displacementIncrement = min(displacementClamp, max(-displacementClamp, displacementIncrement))
                        dofs[dof] += displacementIncrement
                    }
                }

                applyPrescribedValues(&dofs, prescribed: prescribed)

                let postEvaluation = try evaluateGlobal(
                    dofs: dofs,
                    previousStates: states,
                    densities: densities,
                    minimumDensity: clampedMinimumDensity
                )

                if !postEvaluation.residual.allSatisfy(\.isFinite) {
                    throw FEMError.linearSolveFailure("2D explicit solver diverged to non-finite residual values.")
                }

                states = postEvaluation.trialStates
                elementVonMises = postEvaluation.elementVonMises
                elementEnergies = postEvaluation.elementEnergies
                finalReactions = postEvaluation.residual
            }

            let residualNorm = l2Norm(finalReactions, indices: freeDOFs) / max(1.0, sqrt(Float(freeDOFs.count)))
            let loadFactor = Float(step) / Float(problem.controls.loadSteps)
            let maxEquivalentPlastic = states.map(\.equivalentPlasticStrain).max() ?? 0
            let maxDamage = states.map(\.damage).max() ?? 0

            stepHistory.append(
                StepResult2D(
                    step: step,
                    loadFactor: loadFactor,
                    iterations: substeps,
                    residualNorm: residualNorm,
                    maxEquivalentPlasticStrain: maxEquivalentPlastic,
                    maxDamage: maxDamage
                )
            )

            stepSnapshots.append(
                StepSnapshot2D(
                    step: step,
                    loadFactor: loadFactor,
                    displacements: dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count),
                    elementStates: states,
                    elementVonMises: elementVonMises,
                    elementEnergies: elementEnergies,
                    elementDensities: densities
                )
            )
        }

        return SolveResult2D(
            displacements: dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count),
            reactions: finalReactions,
            elementStates: states,
            elementVonMises: elementVonMises,
            elementEnergies: elementEnergies,
            elementDensities: densities,
            stepHistory: stepHistory,
            stepSnapshots: stepSnapshots,
            backendName: backendName,
            converged: true
        )
    }

    private func evaluateGlobal(
        dofs: [Float],
        previousStates: [ElementState],
        densities: [Float],
        minimumDensity: Float
    ) throws -> ExplicitGlobalEvaluation2D {
        let nodalDisplacements = dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count)
        let elementEvaluation = try evaluator.evaluate(
            displacements: nodalDisplacements,
            previousStates: previousStates,
            densities: densities,
            densityPenalty: explicitControls.densityPenalty,
            minimumDensity: minimumDensity
        )

        var residual = Array(repeating: Float(0), count: dofCount)

        for (elementIndex, element) in preparedMesh.elements.enumerated() {
            let forceBase = elementIndex * 3
            for localIndex in 0..<3 {
                let node = Int(element.nodeIDs[localIndex])
                let force = elementEvaluation.elementNodeForces[forceBase + localIndex]
                let dofBase = node * 2
                residual[dofBase] += force.x
                residual[dofBase + 1] += force.y
            }
        }

        return ExplicitGlobalEvaluation2D(
            residual: residual,
            trialStates: elementEvaluation.trialStates,
            elementVonMises: elementEvaluation.elementVonMises,
            elementEnergies: elementEvaluation.elementEnergies
        )
    }

    private func buildLumpedMasses(densities: [Float], minimumDensity: Float) -> [Float] {
        var masses = Array(repeating: Float(0), count: dofCount)

        for (elementIndex, element) in preparedMesh.elements.enumerated() {
            let density = max(minimumDensity, min(1.0, densities[elementIndex]))
            let elementMass = explicitControls.massDensity * density * problem.thickness * element.area
            let nodalMass = elementMass / 3.0

            for localIndex in 0..<3 {
                let node = Int(element.nodeIDs[localIndex])
                let dofBase = node * 2
                masses[dofBase] += nodalMass
                masses[dofBase + 1] += nodalMass
            }
        }

        for index in masses.indices {
            masses[index] = max(masses[index], 1e-8)
        }

        return masses
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
        prescribedDisplacements: [PrescribedDisplacement2D],
        dofCount: Int
    ) throws -> [Int: Float] {
        var map: [Int: Float] = [:]
        for constraint in prescribedDisplacements {
            if !(0..<2).contains(constraint.component) {
                throw FEMError.invalidBoundaryCondition(
                    "2D node \(constraint.node) uses invalid displacement component \(constraint.component)."
                )
            }

            let dof = constraint.dof
            if dof < 0 || dof >= dofCount {
                throw FEMError.invalidBoundaryCondition("2D DOF \(dof) is out of range.")
            }

            if let existing = map[dof], abs(existing - constraint.value) > 1e-8 {
                throw FEMError.invalidBoundaryCondition("Conflicting prescribed displacement values for 2D DOF \(dof).")
            }

            map[dof] = constraint.value
        }
        return map
    }

    private static func makeMetalEvaluator(
        preparedMesh: PreparedMesh2D,
        material: MaterialParameters,
        thickness: Float
    ) throws -> ElementEvaluator2D? {
        #if canImport(Metal)
        return try MetalElementEvaluator2D(preparedMesh: preparedMesh, material: material, thickness: thickness)
        #else
        _ = preparedMesh
        _ = material
        _ = thickness
        return nil
        #endif
    }
}
