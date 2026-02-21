import Foundation

public enum ExampleProblems2D {
    public static func displacementControlledTension(
        length: Float = 1.0,
        height: Float = 0.25,
        nx: Int = 8,
        ny: Int = 3,
        order: ElementOrder2D = .linear,
        subdivisionLevels: Int = 0,
        endDisplacement: Float = 0.08,
        loadSteps: Int = 12
    ) -> FEMProblem2D {
        var mesh = Mesh2D.rectangularPlate(
            length: length,
            height: height,
            nx: nx,
            ny: ny,
            order: order
        )

        if subdivisionLevels > 0 {
            mesh = mesh.subdivided(levels: subdivisionLevels)
        }

        var prescribed: [PrescribedDisplacement2D] = []
        prescribed.reserveCapacity(mesh.nodes.count * 2)

        for (index, node) in mesh.nodes.enumerated() {
            if abs(node.x) < 1e-6 {
                prescribed.append(PrescribedDisplacement2D(node: index, component: 0, value: 0))
                prescribed.append(PrescribedDisplacement2D(node: index, component: 1, value: 0))
            } else if abs(node.x - length) < 1e-6 {
                prescribed.append(PrescribedDisplacement2D(node: index, component: 0, value: endDisplacement))
            }
        }

        let material = MaterialParameters(
            youngsModulus: 70_000,
            poissonRatio: 0.33,
            yieldStress: 250,
            hardeningModulus: 600,
            damageOnset: 0.02,
            damageSlope: 0.25,
            damageCap: 0.95
        )

        let controls = SolverControls(
            loadSteps: max(2, loadSteps),
            maxNewtonIterations: 35,
            residualTolerance: 2e-3,
            finiteDifferenceStep: 2e-4,
            lineSearchFloor: 1.0 / 128.0,
            linearSolverTolerance: 2e-5,
            linearSolverMaxIterations: 900
        )

        return FEMProblem2D(
            mesh: mesh,
            material: material,
            thickness: 1.0,
            prescribedDisplacements: prescribed,
            controls: controls
        )
    }

    public static func simpleShear(
        size: Float = 1.0,
        nx: Int = 6,
        ny: Int = 6,
        order: ElementOrder2D = .linear,
        subdivisionLevels: Int = 0,
        topDisplacement: Float = 0.1,
        loadSteps: Int = 10
    ) -> FEMProblem2D {
        var mesh = Mesh2D.rectangularPlate(
            length: size,
            height: size,
            nx: nx,
            ny: ny,
            order: order
        )

        if subdivisionLevels > 0 {
            mesh = mesh.subdivided(levels: subdivisionLevels)
        }

        var prescribed: [PrescribedDisplacement2D] = []

        for (index, node) in mesh.nodes.enumerated() {
            if abs(node.y) < 1e-6 {
                prescribed.append(PrescribedDisplacement2D(node: index, component: 0, value: 0))
                prescribed.append(PrescribedDisplacement2D(node: index, component: 1, value: 0))
            } else if abs(node.y - size) < 1e-6 {
                prescribed.append(PrescribedDisplacement2D(node: index, component: 0, value: topDisplacement))
                prescribed.append(PrescribedDisplacement2D(node: index, component: 1, value: 0))
            }
        }

        let material = MaterialParameters(
            youngsModulus: 65_000,
            poissonRatio: 0.30,
            yieldStress: 220,
            hardeningModulus: 500,
            damageOnset: 0.025,
            damageSlope: 0.30,
            damageCap: 0.90
        )

        let controls = SolverControls(
            loadSteps: max(2, loadSteps),
            maxNewtonIterations: 35,
            residualTolerance: 2e-3,
            finiteDifferenceStep: 2e-4,
            lineSearchFloor: 1.0 / 128.0,
            linearSolverTolerance: 2e-5,
            linearSolverMaxIterations: 900
        )

        return FEMProblem2D(
            mesh: mesh,
            material: material,
            thickness: 1.0,
            prescribedDisplacements: prescribed,
            controls: controls
        )
    }
}
