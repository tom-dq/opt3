import Foundation

public struct TopologyLiteratureBenchmarkControls2D {
    public var order: ElementOrder2D
    public var integrationScheme: IntegrationScheme2D
    public var subdivisionLevels: Int
    public var iterations: Int
    public var patchRadius: Int
    public var minimumDensity: Float
    public var maximumDensity: Float
    public var moveLimit: Float
    public var referenceElementStride: Int
    public var maxReferenceEvaluationsPerIteration: Int?
    public var objectiveTolerance: Float
    public var densityChangeTolerance: Float
    public var targetVolumeFractionStart: Float?
    public var targetVolumeFractionEnd: Float?
    public var explicitControls: ExplicitSolverControls2D
    public var convoyWorkers: Int
    public var caseIDs: Set<String>?
    public var liangLoads: [Float]
    public var lBracketNX: Int
    public var lBracketNY: Int
    public var uBracketNX: Int
    public var uBracketNY: Int
    public var cantileverNX: Int
    public var cantileverNY: Int
    public var loadSteps: Int
    public var stressWeight: Float
    public var plasticWeight: Float
    public var damageWeight: Float
    public var densityVarianceWeight: Float

    public init(
        order: ElementOrder2D = .linear,
        integrationScheme: IntegrationScheme2D = .full,
        subdivisionLevels: Int = 0,
        iterations: Int = 4,
        patchRadius: Int = 1,
        minimumDensity: Float = 0.001,
        maximumDensity: Float = 1.0,
        moveLimit: Float = 0.15,
        referenceElementStride: Int = 12,
        maxReferenceEvaluationsPerIteration: Int? = 24,
        objectiveTolerance: Float = 1e-7,
        densityChangeTolerance: Float = 1e-4,
        targetVolumeFractionStart: Float? = nil,
        targetVolumeFractionEnd: Float? = nil,
        explicitControls: ExplicitSolverControls2D = ExplicitSolverControls2D(
            stabilization: ExplicitStabilizationControls2D(
                residualBlend: 0.30,
                velocitySmoothing: 0.18,
                displacementSmoothing: 0.02,
                smoothingPasses: 2
            ),
            substepsPerLoadStep: 28,
            relaxationIterationsPerSubstep: 5,
            timeStep: 2.5e-4,
            damping: 0.12,
            massDensity: 2_500,
            densityPenalty: 3.0,
            velocityClamp: 10.0
        ),
        convoyWorkers: Int = min(4, max(1, ProcessInfo.processInfo.activeProcessorCount / 2)),
        caseIDs: Set<String>? = nil,
        liangLoads: [Float] = [5, 10, 15, 20],
        lBracketNX: Int = 24,
        lBracketNY: Int = 24,
        uBracketNX: Int = 28,
        uBracketNY: Int = 20,
        cantileverNX: Int = 36,
        cantileverNY: Int = 12,
        loadSteps: Int = 10,
        stressWeight: Float = 0.05,
        plasticWeight: Float = 0,
        damageWeight: Float = 0,
        densityVarianceWeight: Float = 12
    ) {
        self.order = order
        self.integrationScheme = integrationScheme
        self.subdivisionLevels = max(0, subdivisionLevels)
        self.iterations = max(1, iterations)
        self.patchRadius = max(0, patchRadius)
        self.minimumDensity = max(0.001, min(1.0, minimumDensity))
        self.maximumDensity = max(self.minimumDensity, min(1.0, maximumDensity))
        self.moveLimit = max(1e-3, moveLimit)
        self.referenceElementStride = max(1, referenceElementStride)
        if let maxReferenceEvaluationsPerIteration {
            self.maxReferenceEvaluationsPerIteration = max(1, maxReferenceEvaluationsPerIteration)
        } else {
            self.maxReferenceEvaluationsPerIteration = nil
        }
        self.objectiveTolerance = max(1e-8, objectiveTolerance)
        self.densityChangeTolerance = max(0, densityChangeTolerance)
        self.targetVolumeFractionStart = targetVolumeFractionStart
        self.targetVolumeFractionEnd = targetVolumeFractionEnd
        self.explicitControls = explicitControls
        self.convoyWorkers = max(1, convoyWorkers)
        self.caseIDs = caseIDs
        self.liangLoads = liangLoads.isEmpty ? [10] : liangLoads.map { max(0.5, $0) }.sorted()
        self.lBracketNX = max(8, lBracketNX)
        self.lBracketNY = max(8, lBracketNY)
        self.uBracketNX = max(8, uBracketNX)
        self.uBracketNY = max(8, uBracketNY)
        self.cantileverNX = max(8, cantileverNX)
        self.cantileverNY = max(4, cantileverNY)
        self.loadSteps = max(2, loadSteps)
        self.stressWeight = max(0, stressWeight)
        self.plasticWeight = max(0, plasticWeight)
        self.damageWeight = max(0, damageWeight)
        self.densityVarianceWeight = max(0, densityVarianceWeight)
    }
}

public struct TopologyLiteratureRun2D {
    public var caseID: String
    public var caseName: String
    public var reference: String
    public var profileName: String
    public var loadMagnitude: Float?
    public var elementCount: Int
    public var runtimeSeconds: Double
    public var objective: Float
    public var compliance: Float
    public var averageDensity: Float
    public var finalResidualNorm: Float
    public var maxVonMises: Float
    public var maxEquivalentPlasticStrain: Float
    public var maxDamage: Float
    public var yieldOnsetLoadFactor: Float
    public var rightEdgeDensity: Float
    public var converged: Bool
    public var iterationsCompleted: Int
    public var problem: FEMProblem2D
    public var topologyResult: TopologyOptimizationResult2D
}

public struct TopologyLiteratureTrendCheck2D {
    public var caseName: String
    public var reference: String
    public var metric: String
    public var observed: String
    public var target: String
    public var passed: Bool
}

public struct TopologyLiteratureBenchmarkSweep2D {
    public var runs: [TopologyLiteratureRun2D]
    public var checks: [TopologyLiteratureTrendCheck2D]

    public var allPassed: Bool {
        checks.allSatisfy(\.passed)
    }
}

public enum TopologyLiteratureBenchmarks2D {
    public static func runConvoy(
        backendChoice: ComputeBackendChoice = .auto,
        controls: TopologyLiteratureBenchmarkControls2D = TopologyLiteratureBenchmarkControls2D()
    ) throws -> TopologyLiteratureBenchmarkSweep2D {
        let profiles = objectiveProfiles(controls: controls)
        let cases = buildCases(controls: controls)
        var jobs: [TopologyBenchmarkJob2D] = []

        for benchmark in cases where isEnabled(caseID: benchmark.id, controls: controls) {
            for profile in benchmark.profileIDs {
                guard profiles[profile] != nil else {
                    continue
                }
                jobs.append(TopologyBenchmarkJob2D(caseDefinition: benchmark, profileID: profile))
            }
        }

        if jobs.isEmpty {
            return TopologyLiteratureBenchmarkSweep2D(runs: [], checks: [])
        }

        let workerCount = workerCountForConvoy(backendChoice: backendChoice, requested: controls.convoyWorkers, jobs: jobs.count)
        let queue = OperationQueue()
        queue.name = "mayorfem.topopt.literature.convoy"
        queue.maxConcurrentOperationCount = workerCount

        let lock = NSLock()
        var runs: [TopologyLiteratureRun2D] = []
        var failures: [String] = []

        for job in jobs {
            queue.addOperation {
                do {
                    guard let profile = profiles[job.profileID] else {
                        throw FEMError.invalidBoundaryCondition("Missing objective profile \(job.profileID).")
                    }
                    let run = try executeJob(
                        job: job,
                        profile: profile,
                        backendChoice: backendChoice,
                        controls: controls
                    )
                    lock.lock()
                    runs.append(run)
                    lock.unlock()
                } catch {
                    lock.lock()
                    failures.append("case=\(job.caseDefinition.id), profile=\(job.profileID): \(error)")
                    lock.unlock()
                }
            }
        }

        queue.waitUntilAllOperationsAreFinished()

        if !failures.isEmpty {
            let message = failures.sorted().joined(separator: " | ")
            throw FEMError.linearSolveFailure("Topology literature convoy failed: \(message)")
        }

        runs.sort { lhs, rhs in
            if lhs.caseID == rhs.caseID {
                if lhs.profileName == rhs.profileName {
                    return (lhs.loadMagnitude ?? 0) < (rhs.loadMagnitude ?? 0)
                }
                return lhs.profileName < rhs.profileName
            }
            return lhs.caseID < rhs.caseID
        }

        let checks = buildTrendChecks(from: runs)
        return TopologyLiteratureBenchmarkSweep2D(runs: runs, checks: checks)
    }

    private static func executeJob(
        job: TopologyBenchmarkJob2D,
        profile: ObjectiveProfile2D,
        backendChoice: ComputeBackendChoice,
        controls: TopologyLiteratureBenchmarkControls2D
    ) throws -> TopologyLiteratureRun2D {
        let problem = job.caseDefinition.problemBuilder(controls)

        let targetStart = controls.targetVolumeFractionStart ?? job.caseDefinition.defaultVolumeFraction
        let targetEnd = controls.targetVolumeFractionEnd ?? job.caseDefinition.defaultVolumeFraction
        let explicitControls = tunedExplicitControls(caseID: job.caseDefinition.id, base: controls.explicitControls)

        let topoptControls = TopologyOptimizationControls2D(
            iterations: controls.iterations,
            patchRadius: controls.patchRadius,
            minimumDensity: controls.minimumDensity,
            maximumDensity: controls.maximumDensity,
            targetVolumeFractionStart: targetStart,
            targetVolumeFractionEnd: targetEnd,
            moveLimit: controls.moveLimit,
            referenceElementStride: controls.referenceElementStride,
            maxReferenceEvaluationsPerIteration: controls.maxReferenceEvaluationsPerIteration,
            objectiveTolerance: controls.objectiveTolerance,
            densityChangeTolerance: controls.densityChangeTolerance,
            explicitControls: explicitControls
        )

        let started = Date()
        let optimizer = try PatchTopologyOptimizer2D(
            problem: problem,
            backendChoice: backendChoice,
            controls: topoptControls,
            objective: profile.objective
        )
        let topopt = try optimizer.optimize()
        let elapsed = Date().timeIntervalSince(started)

        let elementCount = topopt.densities.count
        let compliance = topopt.finalSolve.elementEnergies.reduce(0, +)
        let objective = topopt.history.last?.objective ?? compliance
        let averageDensity = topopt.densities.reduce(0, +) / Float(max(1, elementCount))
        let finalResidualNorm = topopt.finalSolve.stepHistory.last?.residualNorm ?? 0
        let maxVonMises = topopt.finalSolve.elementVonMises.max() ?? 0
        let maxEquivalentPlastic = topopt.finalSolve.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
        let maxDamage = topopt.finalSolve.elementStates.map(\.damage).max() ?? 0
        let yieldThreshold = max(5e-4, 0.15 * maxEquivalentPlastic)
        let yieldOnset = yieldOnsetLoadFactor(stepHistory: topopt.finalSolve.stepHistory, threshold: yieldThreshold)
        let rightEdgeDensity = try TopologyLiteratureProblems2D.rightEdgeDensity(problem: problem, densities: topopt.densities)

        return TopologyLiteratureRun2D(
            caseID: job.caseDefinition.id,
            caseName: job.caseDefinition.name,
            reference: job.caseDefinition.reference,
            profileName: profile.name,
            loadMagnitude: job.caseDefinition.loadMagnitude,
            elementCount: elementCount,
            runtimeSeconds: elapsed,
            objective: objective,
            compliance: compliance,
            averageDensity: averageDensity,
            finalResidualNorm: finalResidualNorm,
            maxVonMises: maxVonMises,
            maxEquivalentPlasticStrain: maxEquivalentPlastic,
            maxDamage: maxDamage,
            yieldOnsetLoadFactor: yieldOnset,
            rightEdgeDensity: rightEdgeDensity,
            converged: topopt.finalSolve.converged,
            iterationsCompleted: topopt.history.count,
            problem: problem,
            topologyResult: topopt
        )
    }

    private static func buildTrendChecks(from runs: [TopologyLiteratureRun2D]) -> [TopologyLiteratureTrendCheck2D] {
        var checks: [TopologyLiteratureTrendCheck2D] = []
        checks.reserveCapacity(16)

        if runs.contains(where: { $0.caseID == "amir_l_bracket_top" }) {
            checks += tradeoffChecks(
                runs: runs,
                caseID: "amir_l_bracket_top",
                caseName: "Amir L-bracket top-right load",
                reference: "Amir et al. 2016, Sec. 5.1",
                expectedYieldRetentionMin: 0.95,
                expectedEdgeDensificationMinPercent: 0.2,
                expectedComplianceDriftMaxPercent: 30
            )
        }

        if runs.contains(where: { $0.caseID == "amir_l_bracket_mid" }) {
            checks += tradeoffChecks(
                runs: runs,
                caseID: "amir_l_bracket_mid",
                caseName: "Amir L-bracket mid-right load",
                reference: "Amir et al. 2016, Sec. 5.1",
                expectedYieldRetentionMin: 0.95,
                expectedEdgeDensificationMinPercent: 0.2,
                expectedComplianceDriftMaxPercent: 30
            )
        }

        if runs.contains(where: { $0.caseID == "amir_u_bracket" }) {
            checks += tradeoffChecks(
                runs: runs,
                caseID: "amir_u_bracket",
                caseName: "Amir U-bracket",
                reference: "Amir et al. 2016, Sec. 5.2",
                expectedYieldRetentionMin: 0.95,
                expectedEdgeDensificationMinPercent: 0.2,
                expectedComplianceDriftMaxPercent: 30
            )
        }

        checks += liangLoadSweepChecks(runs: runs)
        return checks
    }

    private static func tradeoffChecks(
        runs: [TopologyLiteratureRun2D],
        caseID: String,
        caseName: String,
        reference: String,
        expectedYieldRetentionMin: Float,
        expectedEdgeDensificationMinPercent: Float,
        expectedComplianceDriftMaxPercent: Float
    ) -> [TopologyLiteratureTrendCheck2D] {
        guard
            let baseline = runs.first(where: { $0.caseID == caseID && $0.profileName == "compliance_baseline" }),
            let guarded = runs.first(where: { $0.caseID == caseID && $0.profileName == "plasticity_guarded" })
        else {
            return [
                TopologyLiteratureTrendCheck2D(
                    caseName: caseName,
                    reference: reference,
                    metric: "pair availability",
                    observed: "missing baseline or guarded run",
                    target: "both profiles must run",
                    passed: false
                )
            ]
        }

        let threshold = max(5e-4, 0.2 * max(baseline.maxEquivalentPlasticStrain, guarded.maxEquivalentPlasticStrain))
        let baselineYield = yieldOnsetLoadFactor(stepHistory: baseline.topologyResult.finalSolve.stepHistory, threshold: threshold)
        let guardedYield = yieldOnsetLoadFactor(stepHistory: guarded.topologyResult.finalSolve.stepHistory, threshold: threshold)
        let yieldRetention = guardedYield / max(1e-6, baselineYield)

        let compliancePenalty = 100.0 * (guarded.compliance - baseline.compliance) / max(1e-6, baseline.compliance)
        let complianceDrift = abs(compliancePenalty)

        let edgeDensityGain = 100.0 * (guarded.rightEdgeDensity - baseline.rightEdgeDensity)
            / max(1e-6, baseline.rightEdgeDensity)

        let yieldCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "yield-onset retention ratio",
            observed: String(format: "%.3f", yieldRetention),
            target: String(format: ">= %.2f", expectedYieldRetentionMin),
            passed: yieldRetention >= expectedYieldRetentionMin
        )

        let complianceCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "compliance drift percent",
            observed: String(format: "%.2f", complianceDrift),
            target: String(format: "<= %.1f", expectedComplianceDriftMaxPercent),
            passed: complianceDrift <= expectedComplianceDriftMaxPercent
        )

        let edgeDensityCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "loaded-edge densification percent",
            observed: String(format: "%.2f", edgeDensityGain),
            target: String(format: ">= %.1f", expectedEdgeDensificationMinPercent),
            passed: edgeDensityGain >= expectedEdgeDensificationMinPercent
        )

        return [yieldCheck, complianceCheck, edgeDensityCheck]
    }

    private static func liangLoadSweepChecks(runs: [TopologyLiteratureRun2D]) -> [TopologyLiteratureTrendCheck2D] {
        let liangRuns = runs
            .filter { $0.caseID.hasPrefix("liang_cantilever_") && $0.profileName == "plasticity_guarded" }
            .sorted { ($0.loadMagnitude ?? 0) < ($1.loadMagnitude ?? 0) }

        guard liangRuns.count >= 2 else {
            return []
        }

        let densityMonotonic = monotonicIncrease(values: liangRuns.map(\.rightEdgeDensity), tolerance: 3e-3)

        let firstDensity = liangRuns.first?.rightEdgeDensity ?? 0
        let lastDensity = liangRuns.last?.rightEdgeDensity ?? firstDensity
        let rightEdgeGain = 100.0 * (lastDensity - firstDensity) / max(1e-6, firstDensity)

        let reference = "Liang et al. 2023 (Machines), Sec. 3.2"
        let caseName = "Liang cantilever load sweep"

        let convergenceCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "all load cases converged",
            observed: liangRuns.allSatisfy(\.converged) ? "yes" : "no",
            target: "yes",
            passed: liangRuns.allSatisfy(\.converged)
        )

        let densityCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "right-edge density monotonic with load",
            observed: densityMonotonic ? "monotonic" : "not monotonic",
            target: "monotonic increase",
            passed: densityMonotonic
        )

        let gainCheck = TopologyLiteratureTrendCheck2D(
            caseName: caseName,
            reference: reference,
            metric: "right-edge density gain percent",
            observed: String(format: "%.2f", rightEdgeGain),
            target: ">= 0.0",
            passed: rightEdgeGain >= 0.0
        )

        return [convergenceCheck, densityCheck, gainCheck]
    }

    private static func buildCases(
        controls: TopologyLiteratureBenchmarkControls2D
    ) -> [TopologyBenchmarkCaseDefinition2D] {
        var cases: [TopologyBenchmarkCaseDefinition2D] = [
            TopologyBenchmarkCaseDefinition2D(
                id: "amir_l_bracket_top",
                name: "Amir L-bracket top-right load",
                reference: "Amir et al. 2016, Sec. 5.1",
                defaultVolumeFraction: 0.45,
                loadMagnitude: nil,
                profileIDs: ["compliance_baseline", "plasticity_guarded"],
                problemBuilder: { cfg in
                    TopologyLiteratureProblems2D.amirLBracket(
                        loadLocation: .topRight,
                        nx: cfg.lBracketNX,
                        ny: cfg.lBracketNY,
                        order: cfg.order,
                        integrationScheme: cfg.integrationScheme,
                        subdivisionLevels: cfg.subdivisionLevels,
                        endDisplacement: 0.052,
                        loadSteps: cfg.loadSteps
                    )
                }
            ),
            TopologyBenchmarkCaseDefinition2D(
                id: "amir_l_bracket_mid",
                name: "Amir L-bracket mid-right load",
                reference: "Amir et al. 2016, Sec. 5.1",
                defaultVolumeFraction: 0.4825,
                loadMagnitude: nil,
                profileIDs: ["compliance_baseline", "plasticity_guarded"],
                problemBuilder: { cfg in
                    TopologyLiteratureProblems2D.amirLBracket(
                        loadLocation: .midRight,
                        nx: cfg.lBracketNX,
                        ny: cfg.lBracketNY,
                        order: cfg.order,
                        integrationScheme: cfg.integrationScheme,
                        subdivisionLevels: cfg.subdivisionLevels,
                        endDisplacement: 0.045,
                        loadSteps: cfg.loadSteps
                    )
                }
            ),
            TopologyBenchmarkCaseDefinition2D(
                id: "amir_u_bracket",
                name: "Amir U-bracket",
                reference: "Amir et al. 2016, Sec. 5.2",
                defaultVolumeFraction: 0.40,
                loadMagnitude: nil,
                profileIDs: ["compliance_baseline", "plasticity_guarded"],
                problemBuilder: { cfg in
                    TopologyLiteratureProblems2D.amirUBracket(
                        nx: cfg.uBracketNX,
                        ny: cfg.uBracketNY,
                        order: cfg.order,
                        integrationScheme: cfg.integrationScheme,
                        subdivisionLevels: cfg.subdivisionLevels,
                        endDisplacement: 0.042,
                        loadSteps: cfg.loadSteps
                    )
                }
            ),
        ]

        for load in controls.liangLoads {
            let label = String(format: "liang_cantilever_%04.1f", load).replacingOccurrences(of: ".", with: "p")
            cases.append(
                TopologyBenchmarkCaseDefinition2D(
                    id: label,
                    name: String(format: "Liang cantilever load %.1f", load),
                    reference: "Liang et al. 2023 (Machines), Sec. 3.2",
                    defaultVolumeFraction: 0.50,
                    loadMagnitude: load,
                    profileIDs: ["plasticity_guarded"],
                    problemBuilder: { cfg in
                        TopologyLiteratureProblems2D.liangCantilever(
                            loadMagnitude: load,
                            nx: cfg.cantileverNX,
                            ny: cfg.cantileverNY,
                            order: cfg.order,
                            integrationScheme: cfg.integrationScheme,
                            subdivisionLevels: cfg.subdivisionLevels,
                            loadSteps: cfg.loadSteps
                        )
                    }
                )
            )
        }

        return cases
    }

    private static func objectiveProfiles(
        controls: TopologyLiteratureBenchmarkControls2D
    ) -> [String: ObjectiveProfile2D] {
        let guardedObjective = plasticityGuardedObjective(
            stressWeight: controls.stressWeight,
            plasticWeight: controls.plasticWeight,
            damageWeight: controls.damageWeight,
            densityVarianceWeight: controls.densityVarianceWeight
        )

        return [
            "compliance_baseline": ObjectiveProfile2D(
                name: "compliance_baseline",
                objective: TopologyObjectives2D.compliance
            ),
            "plasticity_guarded": ObjectiveProfile2D(
                name: "plasticity_guarded",
                objective: guardedObjective
            ),
        ]
    }

    private static func plasticityGuardedObjective(
        stressWeight: Float,
        plasticWeight: Float,
        damageWeight: Float,
        densityVarianceWeight: Float
    ) -> PatchObjectiveFunction2D {
        let safeStressWeight = max(0, stressWeight)
        let safePlasticWeight = max(0, plasticWeight)
        let safeDamageWeight = max(0, damageWeight)
        let safeDensityVarianceWeight = max(0, densityVarianceWeight)

        return { input in
            let compliance = TopologyObjectives2D.compliance(input)
            let meanVonMises = TopologyObjectives2D.meanVonMises(input)

            var meanEquivalentPlastic: Float = 0
            var meanDamage: Float = 0
            var count: Float = 0

            for elementIndex in input.patchElements where elementIndex >= 0 && elementIndex < input.solveResult.elementStates.count {
                let state = input.solveResult.elementStates[elementIndex]
                meanEquivalentPlastic += state.equivalentPlasticStrain
                meanDamage += state.damage
                count += 1
            }

            if count > 0 {
                meanEquivalentPlastic /= count
                meanDamage /= count
            }

            var meanDensity: Float = 0
            var densityCount: Float = 0
            for elementIndex in input.patchElements where elementIndex >= 0 && elementIndex < input.densities.count {
                meanDensity += input.densities[elementIndex]
                densityCount += 1
            }
            if densityCount > 0 {
                meanDensity /= densityCount
            }

            var densityVariance: Float = 0
            if densityCount > 0 {
                for elementIndex in input.patchElements where elementIndex >= 0 && elementIndex < input.densities.count {
                    let centered = input.densities[elementIndex] - meanDensity
                    densityVariance += centered * centered
                }
                densityVariance /= densityCount
            }

            return compliance
                + safeStressWeight * meanVonMises
                + safePlasticWeight * meanEquivalentPlastic
                + safeDamageWeight * meanDamage
                + safeDensityVarianceWeight * densityVariance
        }
    }

    private static func isEnabled(caseID: String, controls: TopologyLiteratureBenchmarkControls2D) -> Bool {
        guard let requestedIDs = controls.caseIDs, !requestedIDs.isEmpty else {
            return true
        }
        if requestedIDs.contains(caseID) {
            return true
        }
        return requestedIDs.contains { prefix in
            caseID.hasPrefix(prefix)
        }
    }

    private static func workerCountForConvoy(
        backendChoice: ComputeBackendChoice,
        requested: Int,
        jobs: Int
    ) -> Int {
        if backendChoice == .metal {
            return 1
        }
        return max(1, min(requested, jobs))
    }

    private static func tunedExplicitControls(
        caseID: String,
        base: ExplicitSolverControls2D
    ) -> ExplicitSolverControls2D {
        if caseID.hasPrefix("liang_cantilever_") {
            return ExplicitSolverControls2D(
                stabilization: ExplicitStabilizationControls2D(
                    residualBlend: max(base.stabilization.residualBlend, 0.45),
                    velocitySmoothing: max(base.stabilization.velocitySmoothing, 0.25),
                    displacementSmoothing: max(base.stabilization.displacementSmoothing, 0.05),
                    smoothingPasses: max(base.stabilization.smoothingPasses, 2)
                ),
                substepsPerLoadStep: max(base.substepsPerLoadStep, 36),
                relaxationIterationsPerSubstep: max(base.relaxationIterationsPerSubstep, 5),
                timeStep: min(base.timeStep, 6e-5),
                damping: max(base.damping, 0.20),
                massDensity: base.massDensity,
                densityPenalty: base.densityPenalty,
                velocityClamp: min(base.velocityClamp, 2.5)
            )
        }
        return base
    }
}

private struct ObjectiveProfile2D {
    var name: String
    var objective: PatchObjectiveFunction2D
}

private struct TopologyBenchmarkCaseDefinition2D {
    var id: String
    var name: String
    var reference: String
    var defaultVolumeFraction: Float
    var loadMagnitude: Float?
    var profileIDs: [String]
    var problemBuilder: (TopologyLiteratureBenchmarkControls2D) -> FEMProblem2D
}

private struct TopologyBenchmarkJob2D {
    var caseDefinition: TopologyBenchmarkCaseDefinition2D
    var profileID: String
}

private func monotonicIncrease(values: [Float], tolerance: Float) -> Bool {
    if values.count < 2 {
        return true
    }
    for index in 1..<values.count {
        if values[index] + tolerance < values[index - 1] {
            return false
        }
    }
    return true
}

private func yieldOnsetLoadFactor(stepHistory: [StepResult2D], threshold: Float) -> Float {
    for step in stepHistory where step.maxEquivalentPlasticStrain >= threshold {
        return step.loadFactor
    }
    return 1.0
}
