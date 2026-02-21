import Foundation
import simd

public enum ExplicitLiteratureBenchmarks2D {
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
            loadSteps: 6
        )
        let solver = try ExplicitFEMSolver2D(
            problem: problem,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        )
        let result = try solver.solve()
        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let tolerance: Float = 8e-4

        return BenchmarkResult2D(
            name: "Explicit 2D J2 elastic range (Simo/Hughes uniaxial anchor)",
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
            loadSteps: 8
        )
        let solver = try ExplicitFEMSolver2D(
            problem: problem,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        )
        let result = try solver.solve()
        let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let threshold: Float = 0.0018

        return BenchmarkResult2D(
            name: "Explicit 2D J2 post-yield response (Belytschko-style loading path)",
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
            loadSteps: 8
        )
        let refined = ExampleProblems2D.displacementControlledTension(
            nx: 5,
            ny: 2,
            order: .linear,
            subdivisionLevels: 1,
            endDisplacement: 0.03,
            loadSteps: 8
        )

        let coarseSolver = try ExplicitFEMSolver2D(
            problem: coarse,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        )
        let refinedSolver = try ExplicitFEMSolver2D(
            problem: refined,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        )

        let coarseResult = try coarseSolver.solve()
        let refinedResult = try refinedSolver.solve()
        let coarseReaction = fixedFaceReaction(problem: coarse, reactions: coarseResult.reactions)
        let refinedReaction = fixedFaceReaction(problem: refined, reactions: refinedResult.reactions)
        let relativeDifference = abs(coarseReaction - refinedReaction) / max(1.0, abs(refinedReaction))
        let tolerance: Float = 0.25

        return BenchmarkResult2D(
            name: "Explicit 2D subdivision consistency (mesh refinement trend)",
            passed: relativeDifference <= tolerance,
            metric: relativeDifference,
            tolerance: tolerance,
            detail: String(format: "reaction_rel_diff=%.6e (expected <= %.6e)", relativeDifference, tolerance)
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
            loadSteps: 8
        )

        var shiftedMesh = base.mesh
        for index in shiftedMesh.nodes.indices {
            shiftedMesh.nodes[index] += SIMD2<Float>(5.75, -2.25)
        }
        let shifted = FEMProblem2D(
            mesh: shiftedMesh,
            material: base.material,
            thickness: base.thickness,
            prescribedDisplacements: base.prescribedDisplacements,
            controls: base.controls
        )

        let baseResult = try ExplicitFEMSolver2D(
            problem: base,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        ).solve()
        let shiftedResult = try ExplicitFEMSolver2D(
            problem: shifted,
            explicitControls: benchmarkControls(),
            backendChoice: backendChoice
        ).solve()

        let baseReaction = fixedFaceReaction(problem: base, reactions: baseResult.reactions)
        let shiftedReaction = fixedFaceReaction(problem: shifted, reactions: shiftedResult.reactions)
        let relativeDifference = abs(baseReaction - shiftedReaction) / max(1.0, abs(baseReaction))
        let tolerance: Float = 0.04

        return BenchmarkResult2D(
            name: "Explicit 2D translational invariance (Bonet/Wood patch anchor)",
            passed: relativeDifference <= tolerance,
            metric: relativeDifference,
            tolerance: tolerance,
            detail: String(format: "reaction_rel_diff=%.6e (expected <= %.6e)", relativeDifference, tolerance)
        )
    }

    public static func benchmarkControls() -> ExplicitSolverControls2D {
        ExplicitSolverControls2D(
            substepsPerLoadStep: 96,
            relaxationIterationsPerSubstep: 8,
            timeStep: 1.0e-4,
            damping: 0.22,
            massDensity: 8_000,
            densityPenalty: 1.0,
            velocityClamp: 2.0
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
