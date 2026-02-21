import Foundation

public enum ExampleProblems {
    public static func displacementControlledTension(
        endDisplacement: Float = 0.08,
        loadSteps: Int = 10
    ) -> FEMProblem {
        let mesh = Mesh.fiveTetraBlock(length: 1.0, width: 0.2, height: 0.2)

        var prescribed: [PrescribedDisplacement] = []
        prescribed.reserveCapacity(mesh.nodes.count * 2)

        for (index, node) in mesh.nodes.enumerated() {
            if abs(node.x) < 1e-6 {
                prescribed.append(PrescribedDisplacement(node: index, component: 0, value: 0))
                prescribed.append(PrescribedDisplacement(node: index, component: 1, value: 0))
                prescribed.append(PrescribedDisplacement(node: index, component: 2, value: 0))
            } else if abs(node.x - 1.0) < 1e-6 {
                prescribed.append(PrescribedDisplacement(node: index, component: 0, value: endDisplacement))
            }
        }

        let material = MaterialParameters(
            youngsModulus: 70_000,
            poissonRatio: 0.33,
            yieldStress: 250,
            hardeningModulus: 600,
            damageOnset: 0.025,
            damageSlope: 0.25,
            damageCap: 0.95
        )

        let controls = SolverControls(
            loadSteps: max(2, loadSteps),
            maxNewtonIterations: 30,
            residualTolerance: 1e-3,
            finiteDifferenceStep: 2e-4,
            lineSearchFloor: 1.0 / 128.0,
            linearSolverTolerance: 1e-5,
            linearSolverMaxIterations: 600
        )

        return FEMProblem(
            mesh: mesh,
            material: material,
            prescribedDisplacements: prescribed,
            controls: controls
        )
    }
}
