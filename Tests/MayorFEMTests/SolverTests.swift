import Foundation
import Testing
@testable import MayorFEM

@Test
func linearAndQuadraticMeshConstructionWorks() throws {
    let linear = Mesh2D.rectangularPlate(nx: 4, ny: 2, order: .linear)
    let quadratic = Mesh2D.rectangularPlate(nx: 4, ny: 2, order: .quadratic)

    #expect(linear.elements.count == 16)
    #expect(quadratic.elements.count == 16)
    #expect(quadratic.nodes.count > linear.nodes.count)

    #expect(linear.elements.allSatisfy {
        if case .tri3 = $0 { return true }
        return false
    })

    #expect(quadratic.elements.allSatisfy {
        if case .tri6 = $0 { return true }
        return false
    })
}

@Test
func subdivisionRefinesTriangleCountByFourPerLevel() throws {
    let base = Mesh2D.rectangularPlate(nx: 3, ny: 1, order: .linear)
    let level1 = base.subdivided(levels: 1)
    let level2 = base.subdivided(levels: 2)

    #expect(level1.elements.count == base.elements.count * 4)
    #expect(level2.elements.count == base.elements.count * 16)
    #expect(level2.nodes.count > level1.nodes.count)
}

@Test
func prepared2DMeshHasPositiveAreas() throws {
    let mesh = Mesh2D.rectangularPlate(nx: 5, ny: 2, order: .quadratic).subdivided(levels: 1)
    let prepared = try mesh.prepare()
    #expect(!prepared.elements.isEmpty)
    for element in prepared.elements {
        #expect(element.area > 0)
    }
}

@Test
func explicitDisplacementControlledTensionConvergesLinearAndQuadratic() throws {
    let controls = ExplicitSolverControls2D(
        substepsPerLoadStep: 40,
        timeStep: 2e-4,
        damping: 0.12,
        massDensity: 2_000,
        densityPenalty: 1.0,
        velocityClamp: 8.0
    )

    let linearProblem = ExampleProblems2D.displacementControlledTension(
        nx: 6,
        ny: 2,
        order: .linear,
        endDisplacement: 0.05,
        loadSteps: 5
    )
    let quadraticProblem = ExampleProblems2D.displacementControlledTension(
        nx: 5,
        ny: 2,
        order: .quadratic,
        endDisplacement: 0.05,
        loadSteps: 5
    )

    let linearResult = try ExplicitFEMSolver2D(
        problem: linearProblem,
        explicitControls: controls,
        backendChoice: .cpu
    ).solve()
    let quadraticResult = try ExplicitFEMSolver2D(
        problem: quadraticProblem,
        explicitControls: controls,
        backendChoice: .cpu
    ).solve()

    #expect(linearResult.converged)
    #expect(quadraticResult.converged)
    #expect(linearResult.stepSnapshots.count == 5)
    #expect(quadraticResult.stepSnapshots.count == 5)
}

@Test
func gradualLoadRampEnforcesBoundaryValuesExplicit() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 6,
        ny: 2,
        order: .linear,
        endDisplacement: 0.04,
        loadSteps: 5
    )
    let controls = ExplicitSolverControls2D(
        substepsPerLoadStep: 32,
        timeStep: 2e-4,
        damping: 0.12,
        massDensity: 2_000,
        densityPenalty: 1.0,
        velocityClamp: 8.0
    )

    let result = try ExplicitFEMSolver2D(problem: problem, explicitControls: controls, backendChoice: .cpu).solve()
    var previousLoad: Float = -1
    let loadedBCs = problem.prescribedDisplacements.filter { $0.component == 0 && abs($0.value) > 1e-8 }

    for snapshot in result.stepSnapshots {
        #expect(snapshot.loadFactor > previousLoad)
        previousLoad = snapshot.loadFactor

        for bc in loadedBCs {
            let expected = bc.value * snapshot.loadFactor
            let actual = snapshot.displacements[bc.node].x
            #expect(abs(expected - actual) < 1e-6)
        }
    }
}

@Test
func explicitBenchmarksPassOnCPU() throws {
    let results = try ExplicitLiteratureBenchmarks2D.runAll(backendChoice: .cpu)
    for result in results {
        #expect(result.passed)
    }
}

@Test
func explicitMetalBenchmarksRunWhenAvailable() throws {
    do {
        let results = try ExplicitLiteratureBenchmarks2D.runAll(backendChoice: .metal)
        for result in results {
            #expect(result.passed)
        }
    } catch FEMError.backendUnavailable {
        // Metal availability is environment-dependent; this still validates behavior.
    }
}

@Test
func cpuMetalConsistencyBenchmarkPassesOrSkips() throws {
    let result = try ExplicitLiteratureBenchmarks2D.cpuMetalConsistencyCheck()
    #expect(result.passed)
}

@Test
func topologyOptimizationChangesDensities() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 3,
        ny: 1,
        order: .linear,
        endDisplacement: 0.03,
        loadSteps: 4
    )
    let controls = TopologyOptimizationControls2D(
        iterations: 2,
        patchRadius: 1,
        minimumDensity: 0.1,
        maximumDensity: 1.0,
        moveLimit: 0.2,
        referenceElementStride: 2,
        explicitControls: ExplicitSolverControls2D(
            substepsPerLoadStep: 16,
            timeStep: 2e-4,
            damping: 0.12,
            massDensity: 2_000,
            densityPenalty: 3.0,
            velocityClamp: 8.0
        )
    )
    let optimizer = try PatchTopologyOptimizer2D(
        problem: problem,
        backendChoice: .cpu,
        controls: controls,
        objective: TopologyObjectives2D.compliance
    )
    let result = try optimizer.optimize()
    let prepared = try problem.mesh.prepare()

    #expect(result.history.count <= 2)
    #expect(result.densities.count == prepared.elements.count)
    #expect(result.history.allSatisfy { $0.objective.isFinite })
    #expect(result.densities.contains { $0 < 0.999 })
    #expect(result.finalSolve.elementDensities.count == result.densities.count)
    #expect(result.densityHistory.count == result.history.count + 1)
}

@Test
func topologyOptimizationRespectsTargetVolumeFraction() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 4,
        ny: 2,
        order: .linear,
        endDisplacement: 0.03,
        loadSteps: 4
    )
    let targetVolume: Float = 0.42
    let controls = TopologyOptimizationControls2D(
        iterations: 3,
        patchRadius: 1,
        minimumDensity: 0.1,
        maximumDensity: 1.0,
        targetVolumeFraction: targetVolume,
        moveLimit: 0.15,
        referenceElementStride: 2,
        objectiveTolerance: 0,
        densityChangeTolerance: 0,
        explicitControls: ExplicitSolverControls2D(
            substepsPerLoadStep: 16,
            timeStep: 2e-4,
            damping: 0.12,
            massDensity: 2_000,
            densityPenalty: 3.0,
            velocityClamp: 8.0
        )
    )

    let optimizer = try PatchTopologyOptimizer2D(
        problem: problem,
        backendChoice: .cpu,
        controls: controls,
        objective: TopologyObjectives2D.compliance
    )
    let result = try optimizer.optimize()
    let averageDensity = result.densities.reduce(0, +) / Float(result.densities.count)

    #expect(abs(averageDensity - targetVolume) < 2e-3)
    #expect(result.history.allSatisfy { abs($0.volumeViolation) < 2e-3 })
}

@Test
func topologyOptimizationCanStopEarlyOnDensityChangeTolerance() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 3,
        ny: 1,
        order: .linear,
        endDisplacement: 0.03,
        loadSteps: 4
    )
    let controls = TopologyOptimizationControls2D(
        iterations: 12,
        patchRadius: 1,
        minimumDensity: 0.1,
        maximumDensity: 1.0,
        moveLimit: 0.01,
        referenceElementStride: 3,
        objectiveTolerance: 0,
        densityChangeTolerance: 0.03,
        explicitControls: ExplicitSolverControls2D(
            substepsPerLoadStep: 12,
            timeStep: 2e-4,
            damping: 0.12,
            massDensity: 2_000,
            densityPenalty: 3.0,
            velocityClamp: 8.0
        )
    )

    let optimizer = try PatchTopologyOptimizer2D(
        problem: problem,
        backendChoice: .cpu,
        controls: controls,
        objective: TopologyObjectives2D.compliance
    )
    let result = try optimizer.optimize()

    #expect(result.converged)
    #expect(result.history.count < 12)
    #expect(result.convergenceReason.contains("density-change tolerance"))
}

@Test
func topologyObjectiveVariantsAreFinite() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 3,
        ny: 1,
        order: .linear,
        endDisplacement: 0.02,
        loadSteps: 3
    )
    let controls = TopologyOptimizationControls2D(
        iterations: 1,
        patchRadius: 1,
        minimumDensity: 0.1,
        maximumDensity: 1.0,
        moveLimit: 0.1,
        referenceElementStride: 2,
        explicitControls: ExplicitSolverControls2D(
            substepsPerLoadStep: 10,
            timeStep: 2e-4,
            damping: 0.12,
            massDensity: 2_000,
            densityPenalty: 3.0,
            velocityClamp: 8.0
        )
    )

    let objectives: [PatchObjectiveFunction2D] = [
        TopologyObjectives2D.compliance,
        TopologyObjectives2D.meanVonMises,
        TopologyObjectives2D.maxDamage,
        TopologyObjectives2D.compliancePlusMeanVonMises(weight: 1e-3),
    ]

    for objective in objectives {
        let optimizer = try PatchTopologyOptimizer2D(
            problem: problem,
            backendChoice: .cpu,
            controls: controls,
            objective: objective
        )
        let result = try optimizer.optimize()
        #expect(result.history.allSatisfy { $0.objective.isFinite })
    }
}

@Test
func sparseBiCGSTABSolvesReferenceSystem() throws {
    var builder = SparseMatrixBuilder(dimension: 3)
    builder.add(row: 0, column: 0, value: 4)
    builder.add(row: 0, column: 1, value: 1)
    builder.add(row: 1, column: 0, value: 2)
    builder.add(row: 1, column: 1, value: 3)
    builder.add(row: 1, column: 2, value: 1)
    builder.add(row: 2, column: 1, value: 1)
    builder.add(row: 2, column: 2, value: 2)

    let matrix = builder.build()
    let expected: [Float] = [1.0, -2.0, 0.5]
    let rhs = matrix.multiply(expected)

    let solved = try solveSparseBiCGSTAB(matrix: matrix, rhs: rhs, tolerance: 1e-6, maxIterations: 200)
    #expect(abs(solved[0] - expected[0]) < 1e-4)
    #expect(abs(solved[1] - expected[1]) < 1e-4)
    #expect(abs(solved[2] - expected[2]) < 1e-4)
}

@Test
func visualizationOutputsVTKAndPNGSeriesFromExplicitSolve() throws {
    let problem = ExampleProblems2D.displacementControlledTension(
        nx: 4,
        ny: 2,
        order: .linear,
        endDisplacement: 0.03,
        loadSteps: 4
    )
    let controls = ExplicitSolverControls2D(
        substepsPerLoadStep: 24,
        timeStep: 2e-4,
        damping: 0.12,
        massDensity: 2_000,
        densityPenalty: 1.0,
        velocityClamp: 8.0
    )
    let result = try ExplicitFEMSolver2D(problem: problem, explicitControls: controls, backendChoice: .cpu).solve()

    let root = FileManager.default.temporaryDirectory
        .appendingPathComponent("mayorfem-2d-viz-\(UUID().uuidString)", isDirectory: true)
    let vtkDir = root.appendingPathComponent("vtk", isDirectory: true)
    let pngDir = root.appendingPathComponent("png", isDirectory: true)

    let vtkFiles = try FEMVisualization.writeVTKSeries(
        problem: problem,
        result: result,
        outputDirectory: vtkDir.path,
        deformationScale: 8
    )

    let pngFiles = try FEMVisualization.writePNGSeries(
        problem: problem,
        result: result,
        outputDirectory: pngDir.path,
        deformationScale: 8,
        imageWidth: 640,
        imageHeight: 420
    )

    #expect(vtkFiles.count == 4)
    #expect(pngFiles.count == 4)
    #expect(FileManager.default.fileExists(atPath: vtkDir.appendingPathComponent("series.pvd").path))

    let firstVTK = try String(contentsOfFile: vtkFiles[0], encoding: .utf8)
    #expect(firstVTK.contains("SCALARS von_mises float 1"))
    #expect(firstVTK.contains("SCALARS density float 1"))

    let firstPNG = try Data(contentsOf: URL(fileURLWithPath: pngFiles[0]))
    #expect(firstPNG.count > 8)
    #expect(firstPNG[0] == 0x89)
    #expect(firstPNG[1] == 0x50)
    #expect(firstPNG[2] == 0x4E)
    #expect(firstPNG[3] == 0x47)

    try? FileManager.default.removeItem(at: root)
}
