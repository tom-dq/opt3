import Foundation

public struct BenchmarkResult2D {
    public var name: String
    public var passed: Bool
    public var metric: Float
    public var tolerance: Float
    public var detail: String
}

public enum LiteratureBenchmarks2D {
    public static func runAll(backendChoice: ComputeBackendChoice = .cpu) throws -> [BenchmarkResult2D] {
        [
            try belowYieldElasticityCheck(backendChoice: backendChoice),
            try aboveYieldPlasticityCheck(backendChoice: backendChoice),
            try subdivisionConsistencyCheck(backendChoice: backendChoice),
            try translationInvarianceCheck(backendChoice: backendChoice),
        ]
    }

    public static func belowYieldElasticityCheck(
        backendChoice: ComputeBackendChoice = .cpu
    ) throws -> BenchmarkResult2D {
        let problem = ExampleProblems2D.displacementControlledTension(
            nx: 6,
            ny: 2,
            order: .linear,
            endDisplacement: 0.0018,
            loadSteps: 5
        )
        let solver = try NonlinearFEMSolver2D(problem: problem, backendChoice: backendChoice)
        let result = try solver.solve()

        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let tolerance: Float = 7e-4

        return BenchmarkResult2D(
            name: "2D J2 elastic range check",
            passed: maxEquivalentPlastic <= tolerance,
            metric: maxEquivalentPlastic,
            tolerance: tolerance,
            detail: String(format: "max_eqp=%.6f (expected <= %.6f)", maxEquivalentPlastic, tolerance)
        )
    }

    public static func aboveYieldPlasticityCheck(
        backendChoice: ComputeBackendChoice = .cpu
    ) throws -> BenchmarkResult2D {
        let problem = ExampleProblems2D.displacementControlledTension(
            nx: 6,
            ny: 2,
            order: .linear,
            endDisplacement: 0.018,
            loadSteps: 6
        )
        let solver = try NonlinearFEMSolver2D(problem: problem, backendChoice: backendChoice)
        let result = try solver.solve()

        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let threshold: Float = 0.002

        return BenchmarkResult2D(
            name: "2D J2 post-yield response check",
            passed: maxEquivalentPlastic >= threshold,
            metric: maxEquivalentPlastic,
            tolerance: threshold,
            detail: String(format: "max_eqp=%.6f (expected >= %.6f)", maxEquivalentPlastic, threshold)
        )
    }

    public static func subdivisionConsistencyCheck(
        backendChoice: ComputeBackendChoice = .cpu
    ) throws -> BenchmarkResult2D {
        let coarse = ExampleProblems2D.displacementControlledTension(
            nx: 5,
            ny: 2,
            order: .linear,
            subdivisionLevels: 0,
            endDisplacement: 0.03,
            loadSteps: 6
        )
        let refined = ExampleProblems2D.displacementControlledTension(
            nx: 5,
            ny: 2,
            order: .linear,
            subdivisionLevels: 1,
            endDisplacement: 0.03,
            loadSteps: 6
        )

        let coarseResult = try NonlinearFEMSolver2D(problem: coarse, backendChoice: backendChoice).solve()
        let refinedResult = try NonlinearFEMSolver2D(problem: refined, backendChoice: backendChoice).solve()

        let coarseReaction = fixedFaceReaction(problem: coarse, reactions: coarseResult.reactions)
        let refinedReaction = fixedFaceReaction(problem: refined, reactions: refinedResult.reactions)

        let relDiff = abs(coarseReaction - refinedReaction) / max(1.0, abs(refinedReaction))
        let tolerance: Float = 0.20

        return BenchmarkResult2D(
            name: "2D subdivision consistency check",
            passed: relDiff <= tolerance,
            metric: relDiff,
            tolerance: tolerance,
            detail: String(format: "reaction_rel_diff=%.6e (expected <= %.6e)", relDiff, tolerance)
        )
    }

    public static func translationInvarianceCheck(
        backendChoice: ComputeBackendChoice = .cpu
    ) throws -> BenchmarkResult2D {
        let base = ExampleProblems2D.displacementControlledTension(
            nx: 6,
            ny: 2,
            order: .linear,
            endDisplacement: 0.03,
            loadSteps: 6
        )
        var translatedMesh = base.mesh
        for index in translatedMesh.nodes.indices {
            translatedMesh.nodes[index] += SIMD2<Float>(8.25, -3.5)
        }

        let translated = FEMProblem2D(
            mesh: translatedMesh,
            material: base.material,
            thickness: base.thickness,
            prescribedDisplacements: base.prescribedDisplacements,
            controls: base.controls
        )

        let baseResult = try NonlinearFEMSolver2D(problem: base, backendChoice: backendChoice).solve()
        let translatedResult = try NonlinearFEMSolver2D(problem: translated, backendChoice: backendChoice).solve()

        let baseReaction = fixedFaceReaction(problem: base, reactions: baseResult.reactions)
        let translatedReaction = fixedFaceReaction(problem: translated, reactions: translatedResult.reactions)

        let relDiff = abs(baseReaction - translatedReaction) / max(1.0, abs(baseReaction))
        let tolerance: Float = 8e-3

        return BenchmarkResult2D(
            name: "2D translational invariance check",
            passed: relDiff <= tolerance,
            metric: relDiff,
            tolerance: tolerance,
            detail: String(format: "reaction_rel_diff=%.6e (expected <= %.6e)", relDiff, tolerance)
        )
    }

    private static func fixedFaceReaction(problem: FEMProblem2D, reactions: [Float]) -> Float {
        var total: Float = 0
        for bc in problem.prescribedDisplacements where bc.component == 0 {
            if abs(bc.value) < 1e-8 {
                total += reactions[bc.dof]
            }
        }
        return total
    }
}
