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
