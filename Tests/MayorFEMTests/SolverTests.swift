import Testing
@testable import MayorFEM
import Foundation

@Test
func preparedMeshVolumesStayPositive() throws {
    let mesh = Mesh.fiveTetraBlock()
    let prepared = try mesh.prepare()
    #expect(prepared.elements.count == 5)
    for element in prepared.elements {
        #expect(element.volume > 0)
    }
}

@Test
func displacementControlledTensionConvergesOnCPU() throws {
    let problem = ExampleProblems.displacementControlledTension(endDisplacement: 0.08, loadSteps: 6)
    let solver = try NonlinearFEMSolver(problem: problem, backendChoice: .cpu)
    let result = try solver.solve()

    #expect(result.converged)
    #expect(result.stepHistory.count == 6)
    #expect(result.stepSnapshots.count == 6)

    let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
    #expect(maxEquivalentPlastic > 0)

    let maxDamage = result.elementStates.map(\.damage).max() ?? 0
    #expect(maxDamage >= 0)
}

@Test
func loadRampAndEnforcedDisplacementsAreAppliedGradually() throws {
    let problem = ExampleProblems.displacementControlledTension(endDisplacement: 0.04, loadSteps: 5)
    let solver = try NonlinearFEMSolver(problem: problem, backendChoice: .cpu)
    let result = try solver.solve()

    #expect(result.stepSnapshots.count == 5)

    let loadedBCs = problem.prescribedDisplacements.filter { $0.component == 0 && abs($0.value) > 1e-8 }
    #expect(!loadedBCs.isEmpty)

    var previousLoad: Float = -1
    for snapshot in result.stepSnapshots {
        #expect(snapshot.loadFactor > previousLoad)
        previousLoad = snapshot.loadFactor

        for bc in loadedBCs {
            let expectedValue = bc.value * snapshot.loadFactor
            let actualValue = snapshot.displacements[bc.node].x
            #expect(abs(actualValue - expectedValue) < 1e-5)
        }
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
func literatureBenchmarksPassOnCPU() throws {
    let results = try LiteratureBenchmarks.runAll(backendChoice: .cpu)
    for result in results {
        #expect(result.passed)
    }
}

@Test
func literatureBenchmarksPassOnMetalWhenAvailable() throws {
    do {
        let results = try LiteratureBenchmarks.runAll(backendChoice: .metal)
        for result in results {
            #expect(result.passed)
        }
    } catch FEMError.backendUnavailable {
        // Metal device not available in this environment.
        return
    }
}

@Test
func vtkVisualizationSeriesContainsStressDeformationAndBCFields() throws {
    let problem = ExampleProblems.displacementControlledTension(endDisplacement: 0.02, loadSteps: 3)
    let solver = try NonlinearFEMSolver(problem: problem, backendChoice: .cpu)
    let result = try solver.solve()

    let outputURL = FileManager.default.temporaryDirectory
        .appendingPathComponent("mayorfem-viz-\(UUID().uuidString)", isDirectory: true)

    let files = try FEMVisualization.writeVTKSeries(
        problem: problem,
        result: result,
        outputDirectory: outputURL.path,
        deformationScale: 12
    )

    #expect(files.count == 3)
    #expect(FileManager.default.fileExists(atPath: outputURL.appendingPathComponent("series.pvd").path))

    let firstContent = try String(contentsOfFile: files[0], encoding: .utf8)
    #expect(firstContent.contains("SCALARS von_mises float 1"))
    #expect(firstContent.contains("VECTORS displacement float"))
    #expect(firstContent.contains("VECTORS prescribed_displacement float"))
    #expect(firstContent.contains("SCALARS enforced_dof_count int 1"))

    try? FileManager.default.removeItem(at: outputURL)
}

@Test
func pngVisualizationSeriesContainsFrames() throws {
    let problem = ExampleProblems.displacementControlledTension(endDisplacement: 0.02, loadSteps: 3)
    let solver = try NonlinearFEMSolver(problem: problem, backendChoice: .cpu)
    let result = try solver.solve()

    let outputURL = FileManager.default.temporaryDirectory
        .appendingPathComponent("mayorfem-png-\(UUID().uuidString)", isDirectory: true)

    let files = try FEMVisualization.writePNGSeries(
        problem: problem,
        result: result,
        outputDirectory: outputURL.path,
        deformationScale: 12,
        imageWidth: 640,
        imageHeight: 480
    )

    #expect(files.count == 3)

    let firstData = try Data(contentsOf: URL(fileURLWithPath: files[0]))
    #expect(firstData.count > 8)
    #expect(firstData[0] == 0x89)
    #expect(firstData[1] == 0x50)
    #expect(firstData[2] == 0x4E)
    #expect(firstData[3] == 0x47)

    try? FileManager.default.removeItem(at: outputURL)
}
