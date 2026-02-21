import Foundation
import simd

private struct GlobalEvaluation2D {
    var residual: [Float]
    var trialStates: [ElementState]
    var elementVonMises: [Float]
    var elementEnergies: [Float]
}

public final class NonlinearFEMSolver2D {
    private let problem: FEMProblem2D
    private let preparedMesh: PreparedMesh2D
    private let linearizer: CPUElementLinearizer2D
    private let backendName: String
    private let dofCount: Int
    private let freeDOFs: [Int]
    private let freeDOFToReducedIndex: [Int]
    private let prescribedFinalValues: [Int: Float]

    public init(problem: FEMProblem2D, backendChoice: ComputeBackendChoice = .auto) throws {
        self.problem = problem
        self.preparedMesh = try problem.mesh.prepare()
        self.dofCount = preparedMesh.nodes.count * 2
        self.linearizer = CPUElementLinearizer2D(
            preparedMesh: preparedMesh,
            material: problem.material,
            thickness: problem.thickness
        )

        let prescribedFinalValues = try NonlinearFEMSolver2D.buildPrescribedMap(
            prescribedDisplacements: problem.prescribedDisplacements,
            dofCount: dofCount
        )
        self.prescribedFinalValues = prescribedFinalValues

        let prescribedSet = Set(prescribedFinalValues.keys)
        self.freeDOFs = (0..<dofCount).filter { !prescribedSet.contains($0) }

        if freeDOFs.isEmpty {
            throw FEMError.noFreeDOFs
        }

        var freeDOFToReducedIndex = Array(repeating: -1, count: dofCount)
        for (reduced, global) in freeDOFs.enumerated() {
            freeDOFToReducedIndex[global] = reduced
        }
        self.freeDOFToReducedIndex = freeDOFToReducedIndex

        switch backendChoice {
        case .metal:
            throw FEMError.backendUnavailable("2D solver Metal backend is not implemented yet.")
        case .cpu, .auto:
            self.backendName = "cpu-2d"
        }
    }

    public func solve() throws -> SolveResult2D {
        var dofs = Array(repeating: Float(0), count: dofCount)
        var states = Array(repeating: ElementState.zero, count: preparedMesh.elements.count)
        var elementVonMises = Array(repeating: Float(0), count: preparedMesh.elements.count)
        var elementEnergies = Array(repeating: Float(0), count: preparedMesh.elements.count)
        var stepHistory: [StepResult2D] = []
        var stepSnapshots: [StepSnapshot2D] = []
        var finalReactions = Array(repeating: Float(0), count: dofCount)

        for step in 1...problem.controls.loadSteps {
            let loadFactor = Float(step) / Float(problem.controls.loadSteps)
            let prescribedAtStep = prescribedValues(for: loadFactor)
            applyPrescribedValues(&dofs, prescribed: prescribedAtStep)

            let baseStates = states
            var convergedEvaluation: GlobalEvaluation2D?
            var iterations = 0
            var residualNorm: Float = .greatestFiniteMagnitude

            for iteration in 1...problem.controls.maxNewtonIterations {
                iterations = iteration
                let evaluation = evaluateGlobal(dofs: dofs, previousStates: baseStates)
                residualNorm = l2Norm(evaluation.residual, indices: freeDOFs)

                if residualNorm <= problem.controls.residualTolerance {
                    convergedEvaluation = evaluation
                    break
                }

                let tangent = assembleSparseTangent(dofs: dofs, previousStates: baseStates)
                let rhs = freeDOFs.map { -evaluation.residual[$0] }
                let increment: [Float]
                do {
                    increment = try solveSparseBiCGSTAB(
                        matrix: tangent,
                        rhs: rhs,
                        tolerance: problem.controls.linearSolverTolerance,
                        maxIterations: problem.controls.linearSolverMaxIterations
                    )
                } catch {
                    increment = try solveDenseLinearSystem(coefficients: tangent.toDenseArray(), rhs: rhs)
                }

                var accepted = false
                var alpha: Float = 1.0
                while alpha >= problem.controls.lineSearchFloor {
                    var candidate = dofs
                    for (index, dofIndex) in freeDOFs.enumerated() {
                        candidate[dofIndex] += alpha * increment[index]
                    }
                    applyPrescribedValues(&candidate, prescribed: prescribedAtStep)

                    let candidateEvaluation = evaluateGlobal(dofs: candidate, previousStates: baseStates)
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
                let terminalEvaluation = evaluateGlobal(dofs: dofs, previousStates: baseStates)
                residualNorm = l2Norm(terminalEvaluation.residual, indices: freeDOFs)
                throw FEMError.newtonFailed(step: step, residualNorm: residualNorm)
            }

            states = convergedEvaluation!.trialStates
            elementVonMises = convergedEvaluation!.elementVonMises
            elementEnergies = convergedEvaluation!.elementEnergies
            finalReactions = convergedEvaluation!.residual

            let maxEquivalentPlastic = states.map(\.equivalentPlasticStrain).max() ?? 0
            let maxDamage = states.map(\.damage).max() ?? 0
            stepHistory.append(
                StepResult2D(
                    step: step,
                    loadFactor: loadFactor,
                    iterations: iterations,
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
                    elementDensities: Array(repeating: 1.0, count: preparedMesh.elements.count)
                )
            )
        }

        return SolveResult2D(
            displacements: dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count),
            reactions: finalReactions,
            elementStates: states,
            elementVonMises: elementVonMises,
            elementEnergies: elementEnergies,
            elementDensities: Array(repeating: 1.0, count: preparedMesh.elements.count),
            stepHistory: stepHistory,
            stepSnapshots: stepSnapshots,
            backendName: backendName,
            converged: true
        )
    }

    private func evaluateGlobal(dofs: [Float], previousStates: [ElementState]) -> GlobalEvaluation2D {
        let nodalDisplacements = dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count)

        var residual = Array(repeating: Float(0), count: dofCount)
        var trialStates = Array(repeating: ElementState.zero, count: preparedMesh.elements.count)
        var vonMises = Array(repeating: Float(0), count: preparedMesh.elements.count)
        var energies = Array(repeating: Float(0), count: preparedMesh.elements.count)

        for elementIndex in preparedMesh.elements.indices {
            let geometry = preparedMesh.elements[elementIndex]
            let reference = gatherElementNodalVectors2D(geometry: geometry, values: preparedMesh.nodes)
            let localDisplacements = gatherElementNodalVectors2D(geometry: geometry, values: nodalDisplacements)

            let response = computeLocalElementResponse2D(
                geometry: geometry,
                reference: reference,
                displacements: localDisplacements,
                previousState: previousStates[elementIndex],
                material: problem.material,
                thickness: problem.thickness
            )

            trialStates[elementIndex] = response.updatedState
            vonMises[elementIndex] = response.vonMisesStress
            energies[elementIndex] = response.strainEnergy

            let nodeIDs = geometry.nodeIDs
            for localIndex in 0..<3 {
                let node = Int(nodeIDs[localIndex])
                let force = response.forces[localIndex]
                let dofBase = node * 2
                residual[dofBase] += force.x
                residual[dofBase + 1] += force.y
            }
        }

        return GlobalEvaluation2D(
            residual: residual,
            trialStates: trialStates,
            elementVonMises: vonMises,
            elementEnergies: energies
        )
    }

    private func assembleSparseTangent(dofs: [Float], previousStates: [ElementState]) -> SparseMatrix {
        let nodalDisplacements = dofVectorToNodalDisplacements2D(dofs, nodeCount: preparedMesh.nodes.count)
        var builder = SparseMatrixBuilder(
            dimension: freeDOFs.count,
            reserve: preparedMesh.elements.count * 6 * 6
        )

        for elementIndex in preparedMesh.elements.indices {
            let linearization = linearizer.linearize(
                elementIndex: elementIndex,
                displacements: nodalDisplacements,
                previousState: previousStates[elementIndex],
                finiteDifferenceStep: problem.controls.finiteDifferenceStep
            )

            let localDOFs = localDOFIndices(for: linearization.nodeIDs)

            for localRow in 0..<6 {
                let reducedRow = freeDOFToReducedIndex[localDOFs[localRow]]
                if reducedRow < 0 {
                    continue
                }

                for localColumn in 0..<6 {
                    let reducedColumn = freeDOFToReducedIndex[localDOFs[localColumn]]
                    if reducedColumn < 0 {
                        continue
                    }

                    let value = linearization.localTangent[localRow * 6 + localColumn]
                    builder.add(row: reducedRow, column: reducedColumn, value: value)
                }
            }
        }

        return builder.build()
    }

    private func localDOFIndices(for nodeIDs: SIMD3<UInt32>) -> [Int] {
        let n0 = Int(nodeIDs[0])
        let n1 = Int(nodeIDs[1])
        let n2 = Int(nodeIDs[2])

        return [
            n0 * 2, n0 * 2 + 1,
            n1 * 2, n1 * 2 + 1,
            n2 * 2, n2 * 2 + 1,
        ]
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
}
