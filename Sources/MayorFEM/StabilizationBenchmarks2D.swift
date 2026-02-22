import Foundation
import simd

public struct StabilizationProfile2D {
    public var name: String
    public var controls: ExplicitStabilizationControls2D
}

public struct StabilizationBenchmarkCase2D {
    public var name: String
    public var coarseProblem: FEMProblem2D
    public var refinedProblem: FEMProblem2D
}

public struct StabilizationBenchmarkRow2D {
    public var caseName: String
    public var profileName: String
    public var roughnessIndex: Float
    public var refinementMismatch: Float
    public var finalResidualNorm: Float
    public var maxDisplacementMagnitude: Float
    public var runtimeSeconds: Double
    public var roughnessReductionPercent: Float
    public var mismatchReductionPercent: Float
    public var combinedEffectivenessPercent: Float
}

public struct StabilizationBenchmarkSweep2D {
    public var rows: [StabilizationBenchmarkRow2D]
}

public enum StabilizationBenchmarks2D {
    public static func defaultProfiles() -> [StabilizationProfile2D] {
        [
            StabilizationProfile2D(
                name: "none",
                controls: ExplicitStabilizationControls2D()
            ),
            StabilizationProfile2D(
                name: "residual_filter",
                controls: ExplicitStabilizationControls2D(
                    residualBlend: 0.55
                )
            ),
            StabilizationProfile2D(
                name: "velocity_smooth",
                controls: ExplicitStabilizationControls2D(
                    velocitySmoothing: 0.22,
                    smoothingPasses: 1
                )
            ),
            StabilizationProfile2D(
                name: "displacement_smooth",
                controls: ExplicitStabilizationControls2D(
                    displacementSmoothing: 0.04,
                    smoothingPasses: 1
                )
            ),
            StabilizationProfile2D(
                name: "hybrid_balanced",
                controls: ExplicitStabilizationControls2D(
                    residualBlend: 0.30,
                    velocitySmoothing: 0.18,
                    displacementSmoothing: 0.0,
                    smoothingPasses: 2
                )
            ),
            StabilizationProfile2D(
                name: "hybrid_strong",
                controls: ExplicitStabilizationControls2D(
                    residualBlend: 0.45,
                    velocitySmoothing: 0.28,
                    displacementSmoothing: 0.04,
                    smoothingPasses: 2
                )
            ),
        ]
    }

    public static func defaultCases() -> [StabilizationBenchmarkCase2D] {
        let tensionReduced = StabilizationBenchmarkCase2D(
            name: "tension_reduced",
            coarseProblem: ExampleProblems2D.displacementControlledTension(
                nx: 12,
                ny: 4,
                order: .linear,
                integrationScheme: .reduced,
                subdivisionLevels: 0,
                endDisplacement: 0.06,
                loadSteps: 8
            ),
            refinedProblem: ExampleProblems2D.displacementControlledTension(
                nx: 12,
                ny: 4,
                order: .linear,
                integrationScheme: .reduced,
                subdivisionLevels: 1,
                endDisplacement: 0.06,
                loadSteps: 8
            )
        )

        let shearReduced = StabilizationBenchmarkCase2D(
            name: "shear_reduced",
            coarseProblem: ExampleProblems2D.simpleShear(
                size: 1.0,
                nx: 10,
                ny: 10,
                order: .linear,
                integrationScheme: .reduced,
                subdivisionLevels: 0,
                topDisplacement: 0.08,
                loadSteps: 8
            ),
            refinedProblem: ExampleProblems2D.simpleShear(
                size: 1.0,
                nx: 10,
                ny: 10,
                order: .linear,
                integrationScheme: .reduced,
                subdivisionLevels: 1,
                topDisplacement: 0.08,
                loadSteps: 8
            )
        )

        let tensionFull = StabilizationBenchmarkCase2D(
            name: "tension_full",
            coarseProblem: ExampleProblems2D.displacementControlledTension(
                nx: 12,
                ny: 4,
                order: .linear,
                integrationScheme: .full,
                subdivisionLevels: 0,
                endDisplacement: 0.06,
                loadSteps: 8
            ),
            refinedProblem: ExampleProblems2D.displacementControlledTension(
                nx: 12,
                ny: 4,
                order: .linear,
                integrationScheme: .full,
                subdivisionLevels: 1,
                endDisplacement: 0.06,
                loadSteps: 8
            )
        )

        return [tensionReduced, shearReduced, tensionFull]
    }

    public static func runSweep(
        backendChoice: ComputeBackendChoice = .cpu,
        profiles: [StabilizationProfile2D] = defaultProfiles(),
        cases: [StabilizationBenchmarkCase2D] = defaultCases()
    ) throws -> StabilizationBenchmarkSweep2D {
        var rows: [StabilizationBenchmarkRow2D] = []
        rows.reserveCapacity(profiles.count * cases.count)

        for benchmarkCase in cases {
            for profile in profiles {
                let controls = benchmarkControls(stabilization: profile.controls)
                let started = Date()

                let coarseResult = try ExplicitFEMSolver2D(
                    problem: benchmarkCase.coarseProblem,
                    explicitControls: controls,
                    backendChoice: backendChoice
                ).solve()
                let refinedResult = try ExplicitFEMSolver2D(
                    problem: benchmarkCase.refinedProblem,
                    explicitControls: controls,
                    backendChoice: backendChoice
                ).solve()

                let elapsed = Date().timeIntervalSince(started)
                let roughness = try displacementRoughness(
                    problem: benchmarkCase.coarseProblem,
                    result: coarseResult
                )
                let mismatch = try refinementMismatch(
                    coarseProblem: benchmarkCase.coarseProblem,
                    coarseResult: coarseResult,
                    refinedProblem: benchmarkCase.refinedProblem,
                    refinedResult: refinedResult
                )
                let finalResidual = coarseResult.stepHistory.last?.residualNorm ?? 0
                let maxDisp = coarseResult.displacements.map(simd_length).max() ?? 0

                rows.append(
                    StabilizationBenchmarkRow2D(
                        caseName: benchmarkCase.name,
                        profileName: profile.name,
                        roughnessIndex: roughness,
                        refinementMismatch: mismatch,
                        finalResidualNorm: finalResidual,
                        maxDisplacementMagnitude: maxDisp,
                        runtimeSeconds: elapsed,
                        roughnessReductionPercent: 0,
                        mismatchReductionPercent: 0,
                        combinedEffectivenessPercent: 0
                    )
                )
            }
        }

        let groupedByCase = Dictionary(grouping: rows.indices, by: { rows[$0].caseName })
        for (_, indices) in groupedByCase {
            guard let baselineIndex = indices.first(where: { rows[$0].profileName == "none" }) else {
                continue
            }
            let baselineRoughness = max(rows[baselineIndex].roughnessIndex, 1e-8)
            let baselineMismatch = max(rows[baselineIndex].refinementMismatch, 1e-8)

            for index in indices {
                let roughnessReduction = 100.0 * (baselineRoughness - rows[index].roughnessIndex) / baselineRoughness
                let mismatchReduction = 100.0 * (baselineMismatch - rows[index].refinementMismatch) / baselineMismatch
                rows[index].roughnessReductionPercent = roughnessReduction
                rows[index].mismatchReductionPercent = mismatchReduction
                rows[index].combinedEffectivenessPercent = 0.5 * (roughnessReduction + mismatchReduction)
            }
        }

        rows.sort {
            if $0.caseName == $1.caseName {
                return $0.combinedEffectivenessPercent > $1.combinedEffectivenessPercent
            }
            return $0.caseName < $1.caseName
        }

        return StabilizationBenchmarkSweep2D(rows: rows)
    }

    public static func benchmarkControls(
        stabilization: ExplicitStabilizationControls2D
    ) -> ExplicitSolverControls2D {
        ExplicitSolverControls2D(
            stabilization: stabilization,
            substepsPerLoadStep: 32,
            relaxationIterationsPerSubstep: 4,
            timeStep: 2.5e-4,
            damping: 0.10,
            massDensity: 2_500,
            densityPenalty: 1.0,
            velocityClamp: 10.0
        )
    }

    private static func displacementRoughness(
        problem: FEMProblem2D,
        result: SolveResult2D
    ) throws -> Float {
        let prepared = try problem.mesh.prepare()
        let adjacency = buildNodeAdjacency(preparedMesh: prepared)
        let displacements = result.displacements

        var sumLaplacianSquared: Float = 0
        var sumDisplacementSquared: Float = 0
        var usedNodes: Float = 0

        for node in 0..<prepared.nodes.count {
            let neighbors = adjacency[node]
            if neighbors.isEmpty {
                continue
            }

            var average = SIMD2<Float>.zero
            for neighbor in neighbors {
                average += displacements[neighbor]
            }
            average /= Float(neighbors.count)

            let delta = displacements[node] - average
            sumLaplacianSquared += simd_dot(delta, delta)
            sumDisplacementSquared += simd_dot(displacements[node], displacements[node])
            usedNodes += 1
        }

        if usedNodes <= 0 {
            return 0
        }

        let laplacianRMS = sqrt(sumLaplacianSquared / usedNodes)
        let displacementRMS = sqrt(sumDisplacementSquared / usedNodes)
        return laplacianRMS / max(1e-6, displacementRMS)
    }

    private static func refinementMismatch(
        coarseProblem: FEMProblem2D,
        coarseResult: SolveResult2D,
        refinedProblem: FEMProblem2D,
        refinedResult: SolveResult2D
    ) throws -> Float {
        let coarseNodes = try coarseProblem.mesh.prepare().nodes
        let refinedNodes = try refinedProblem.mesh.prepare().nodes

        var refinedMap: [String: Int] = [:]
        refinedMap.reserveCapacity(refinedNodes.count)
        for (index, node) in refinedNodes.enumerated() {
            refinedMap[coordinateKey(node)] = index
        }

        var sumDiffSquared: Float = 0
        var sumRefSquared: Float = 0
        var matchedCount: Float = 0

        for (coarseIndex, node) in coarseNodes.enumerated() {
            guard let refinedIndex = refinedMap[coordinateKey(node)] else {
                continue
            }

            let coarseDisp = coarseResult.displacements[coarseIndex]
            let refinedDisp = refinedResult.displacements[refinedIndex]
            let diff = coarseDisp - refinedDisp

            sumDiffSquared += simd_dot(diff, diff)
            sumRefSquared += simd_dot(refinedDisp, refinedDisp)
            matchedCount += 1
        }

        if matchedCount < Float(max(4, coarseNodes.count / 2)) {
            throw FEMError.invalidMesh(
                "Failed to map enough coarse nodes into refined mesh during stabilization refinement comparison."
            )
        }

        let diffRMS = sqrt(sumDiffSquared / matchedCount)
        let refRMS = sqrt(sumRefSquared / matchedCount)
        return diffRMS / max(1e-6, refRMS)
    }

    private static func buildNodeAdjacency(preparedMesh: PreparedMesh2D) -> [[Int]] {
        var adjacency = Array(repeating: Set<Int>(), count: preparedMesh.nodes.count)
        for element in preparedMesh.elements {
            let nodeIDs = [
                Int(element.nodeIDs[0]),
                Int(element.nodeIDs[1]),
                Int(element.nodeIDs[2]),
                Int(element.nodeIDs[3]),
            ]
            for i in 0..<4 {
                let a = nodeIDs[i]
                let b = nodeIDs[(i + 1) % 4]
                adjacency[a].insert(b)
                adjacency[b].insert(a)
            }
        }
        return adjacency.map { $0.sorted() }
    }

    private static func coordinateKey(_ node: SIMD2<Float>) -> String {
        let x = Int64((Double(node.x) * 1_000_000.0).rounded())
        let y = Int64((Double(node.y) * 1_000_000.0).rounded())
        return "\(x):\(y)"
    }
}
