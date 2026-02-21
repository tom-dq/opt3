import Foundation
import simd

public struct BenchmarkResult {
    public var name: String
    public var passed: Bool
    public var metric: Float
    public var tolerance: Float
    public var detail: String
}

public enum LiteratureBenchmarks {
    public static func runAll(backendChoice: ComputeBackendChoice = .cpu) throws -> [BenchmarkResult] {
        [
            try belowYieldElasticityCheck(backendChoice: backendChoice),
            try aboveYieldPlasticityCheck(backendChoice: backendChoice),
            try translationInvarianceCheck(backendChoice: backendChoice),
        ]
    }

    public static func belowYieldElasticityCheck(backendChoice: ComputeBackendChoice = .cpu) throws -> BenchmarkResult {
        let problem = stabilized(
            problem: ExampleProblems.displacementControlledTension(endDisplacement: 0.0018, loadSteps: 4)
        )
        let solver = try NonlinearFEMSolver(problem: problem, backendChoice: backendChoice)
        let result = try solver.solve()

        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let tolerance: Float = 5e-4

        return BenchmarkResult(
            name: "J2 elastic range (Simo/Hughes-style uniaxial check)",
            passed: maxEquivalentPlastic <= tolerance,
            metric: maxEquivalentPlastic,
            tolerance: tolerance,
            detail: String(format: "max_eqp=%.6f (expected <= %.6f)", maxEquivalentPlastic, tolerance)
        )
    }

    public static func aboveYieldPlasticityCheck(backendChoice: ComputeBackendChoice = .cpu) throws -> BenchmarkResult {
        let problem = stabilized(
            problem: ExampleProblems.displacementControlledTension(endDisplacement: 0.012, loadSteps: 6)
        )
        let solver = try NonlinearFEMSolver(problem: problem, backendChoice: backendChoice)
        let result = try solver.solve()

        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let threshold: Float = 0.001

        return BenchmarkResult(
            name: "J2 post-yield response (Simo/Hughes-style uniaxial check)",
            passed: maxEquivalentPlastic >= threshold,
            metric: maxEquivalentPlastic,
            tolerance: threshold,
            detail: String(format: "max_eqp=%.6f (expected >= %.6f)", maxEquivalentPlastic, threshold)
        )
    }

    public static func translationInvarianceCheck(backendChoice: ComputeBackendChoice = .cpu) throws -> BenchmarkResult {
        let baseProblem = stabilized(
            problem: ExampleProblems.displacementControlledTension(endDisplacement: 0.03, loadSteps: 6)
        )
        let translatedProblem = translated(problem: baseProblem, by: SIMD3<Float>(7.5, -4.0, 2.25))

        let baseSolver = try NonlinearFEMSolver(problem: baseProblem, backendChoice: backendChoice)
        let translatedSolver = try NonlinearFEMSolver(problem: translatedProblem, backendChoice: backendChoice)

        let baseResult = try baseSolver.solve()
        let translatedResult = try translatedSolver.solve()

        let baseReaction = supportXReaction(problem: baseProblem, reactions: baseResult.reactions)
        let translatedReaction = supportXReaction(problem: translatedProblem, reactions: translatedResult.reactions)

        let denominator = max(1.0, abs(baseReaction))
        let relativeDifference = abs(baseReaction - translatedReaction) / denominator
        let tolerance: Float = 5e-3

        return BenchmarkResult(
            name: "Finite-strain translational invariance (Bonet/Wood patch check)",
            passed: relativeDifference <= tolerance,
            metric: relativeDifference,
            tolerance: tolerance,
            detail: String(
                format: "reaction_rel_diff=%.6e (expected <= %.6e)",
                relativeDifference,
                tolerance
            )
        )
    }

    private static func translated(problem: FEMProblem, by shift: SIMD3<Float>) -> FEMProblem {
        var shiftedMesh = problem.mesh
        for index in shiftedMesh.nodes.indices {
            shiftedMesh.nodes[index] += shift
        }

        return FEMProblem(
            mesh: shiftedMesh,
            material: problem.material,
            prescribedDisplacements: problem.prescribedDisplacements,
            controls: problem.controls
        )
    }

    private static func supportXReaction(problem: FEMProblem, reactions: [Float]) -> Float {
        var total: Float = 0

        for bc in problem.prescribedDisplacements where bc.component == 0 {
            if abs(bc.value) < 1e-8 {
                total += reactions[bc.dof]
            }
        }

        return total
    }

    private static func stabilized(problem: FEMProblem) -> FEMProblem {
        var controls = problem.controls
        controls.maxNewtonIterations = max(controls.maxNewtonIterations, 40)
        controls.residualTolerance = max(controls.residualTolerance, 5e-3)
        controls.linearSolverTolerance = max(controls.linearSolverTolerance, 5e-5)

        return FEMProblem(
            mesh: problem.mesh,
            material: problem.material,
            prescribedDisplacements: problem.prescribedDisplacements,
            controls: controls
        )
    }
}
