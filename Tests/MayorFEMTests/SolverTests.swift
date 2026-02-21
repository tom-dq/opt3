import Testing
@testable import MayorFEM

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

    let maxEquivalentPlastic = result.elementStates.map(\.equivalentPlasticStrain).max() ?? 0
    #expect(maxEquivalentPlastic > 0)

    let maxDamage = result.elementStates.map(\.damage).max() ?? 0
    #expect(maxDamage >= 0)
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
