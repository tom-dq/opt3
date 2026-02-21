import Foundation
import MayorFEM

private enum SolverMode {
    case implicitNewton
    case explicitDynamics
}

private struct CLIOptions {
    var steps: Int = 12
    var displacement: Float = 0.08
    var backend: ComputeBackendChoice = .auto
    var solverMode: SolverMode = .explicitDynamics
    var runBenchmarks: Bool = false
    var visualizationDirectory: String?
    var deformationScale: Float = 10.0
    var imageWidth: Int = 1400
    var imageHeight: Int = 900
    var nx: Int = 8
    var ny: Int = 3
    var order: ElementOrder2D = .linear
    var subdivisionLevels: Int = 0

    var explicitSubsteps: Int = 48
    var explicitRelaxIterations: Int = 6
    var explicitTimeStep: Float = 4e-4
    var explicitDamping: Float = 0.08
    var explicitMassDensity: Float = 1_250.0
    var explicitDensityPenalty: Float = 3.0
    var explicitVelocityClamp: Float = 25.0

    var runTopologyOptimization: Bool = false
    var topologyIterations: Int = 8
    var topologyPatchRadius: Int = 1
    var topologyMinimumDensity: Float = 0.05
    var topologyMaximumDensity: Float = 1.0
    var topologyTargetVolumeFraction: Float?
    var topologyMoveLimit: Float = 0.12
    var topologyReferenceStride: Int = 1
    var topologyObjectiveTolerance: Float = 1e-5
    var topologyDensityChangeTolerance: Float = 1e-4
    var topologyExportEvery: Int = 0
    var objectiveName: String = "compliance"
    var objectiveWeight: Float = 1e-3

    static func parse(arguments: [String]) throws -> CLIOptions {
        var options = CLIOptions()
        var index = 1

        while index < arguments.count {
            let argument = arguments[index]
            switch argument {
            case "--steps":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 1 else {
                    throw FEMError.invalidBoundaryCondition("--steps requires an integer > 1")
                }
                options.steps = value
            case "--disp":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--disp requires a positive floating point value")
                }
                options.displacement = value
            case "--backend":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--backend requires one of auto|metal|cpu")
                }
                switch arguments[index] {
                case "auto": options.backend = .auto
                case "metal": options.backend = .metal
                case "cpu": options.backend = .cpu
                default:
                    throw FEMError.invalidBoundaryCondition("--backend requires one of auto|metal|cpu")
                }
            case "--solver":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--solver requires one of explicit|implicit")
                }
                switch arguments[index] {
                case "explicit": options.solverMode = .explicitDynamics
                case "implicit": options.solverMode = .implicitNewton
                default:
                    throw FEMError.invalidBoundaryCondition("--solver requires one of explicit|implicit")
                }
            case "--benchmarks":
                options.runBenchmarks = true
            case "--visualize":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--visualize requires an output directory")
                }
                options.visualizationDirectory = arguments[index]
            case "--deformation-scale":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--deformation-scale requires a positive value")
                }
                options.deformationScale = value
            case "--image-width":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 256 else {
                    throw FEMError.invalidBoundaryCondition("--image-width requires integer >= 256")
                }
                options.imageWidth = value
            case "--image-height":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 256 else {
                    throw FEMError.invalidBoundaryCondition("--image-height requires integer >= 256")
                }
                options.imageHeight = value
            case "--nx":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--nx requires integer > 0")
                }
                options.nx = value
            case "--ny":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--ny requires integer > 0")
                }
                options.ny = value
            case "--order":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--order requires linear|quadratic")
                }
                switch arguments[index] {
                case "linear": options.order = .linear
                case "quadratic": options.order = .quadratic
                default:
                    throw FEMError.invalidBoundaryCondition("--order requires linear|quadratic")
                }
            case "--subdivide":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--subdivide requires integer >= 0")
                }
                options.subdivisionLevels = value
            case "--explicit-substeps":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--explicit-substeps requires integer > 0")
                }
                options.explicitSubsteps = value
            case "--relax-iters":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--relax-iters requires integer > 0")
                }
                options.explicitRelaxIterations = value
            case "--dt":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--dt requires a positive value")
                }
                options.explicitTimeStep = value
            case "--damping":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--damping requires a non-negative value")
                }
                options.explicitDamping = value
            case "--mass-density":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--mass-density requires a positive value")
                }
                options.explicitMassDensity = value
            case "--density-penalty":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 1 else {
                    throw FEMError.invalidBoundaryCondition("--density-penalty requires value >= 1")
                }
                options.explicitDensityPenalty = value
            case "--velocity-clamp":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--velocity-clamp requires a positive value")
                }
                options.explicitVelocityClamp = value
            case "--topopt":
                options.runTopologyOptimization = true
            case "--topopt-iters":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--topopt-iters requires integer > 0")
                }
                options.topologyIterations = value
            case "--patch-radius":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--patch-radius requires integer >= 0")
                }
                options.topologyPatchRadius = value
            case "--min-density":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0, value <= 1 else {
                    throw FEMError.invalidBoundaryCondition("--min-density requires 0 < value <= 1")
                }
                options.topologyMinimumDensity = value
            case "--max-density":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0, value <= 1 else {
                    throw FEMError.invalidBoundaryCondition("--max-density requires 0 < value <= 1")
                }
                options.topologyMaximumDensity = value
            case "--volume-fraction":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0, value <= 1 else {
                    throw FEMError.invalidBoundaryCondition("--volume-fraction requires 0 < value <= 1")
                }
                options.topologyTargetVolumeFraction = value
            case "--move-limit":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--move-limit requires value > 0")
                }
                options.topologyMoveLimit = value
            case "--reference-stride":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--reference-stride requires integer > 0")
                }
                options.topologyReferenceStride = value
            case "--objective-tol":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--objective-tol requires value >= 0")
                }
                options.topologyObjectiveTolerance = value
            case "--density-change-tol":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--density-change-tol requires value >= 0")
                }
                options.topologyDensityChangeTolerance = value
            case "--topopt-export-every":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--topopt-export-every requires integer >= 0")
                }
                options.topologyExportEvery = value
            case "--objective":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--objective requires a name")
                }
                options.objectiveName = arguments[index]
            case "--objective-weight":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--objective-weight requires value >= 0")
                }
                options.objectiveWeight = value
            case "--help", "-h":
                printUsage()
                exit(0)
            default:
                throw FEMError.invalidBoundaryCondition("Unknown argument: \(argument)")
            }
            index += 1
        }

        return options
    }

    func explicitControls() -> ExplicitSolverControls2D {
        ExplicitSolverControls2D(
            substepsPerLoadStep: explicitSubsteps,
            relaxationIterationsPerSubstep: explicitRelaxIterations,
            timeStep: explicitTimeStep,
            damping: explicitDamping,
            massDensity: explicitMassDensity,
            densityPenalty: explicitDensityPenalty,
            velocityClamp: explicitVelocityClamp
        )
    }

    func topologyControls() -> TopologyOptimizationControls2D {
        TopologyOptimizationControls2D(
            iterations: topologyIterations,
            patchRadius: topologyPatchRadius,
            minimumDensity: topologyMinimumDensity,
            maximumDensity: topologyMaximumDensity,
            targetVolumeFraction: topologyTargetVolumeFraction,
            moveLimit: topologyMoveLimit,
            referenceElementStride: topologyReferenceStride,
            objectiveTolerance: topologyObjectiveTolerance,
            densityChangeTolerance: topologyDensityChangeTolerance,
            explicitControls: explicitControls()
        )
    }

    static func printUsage() {
        print(
            """
            mayor-fem (2D)

            Usage:
              swift run mayor-fem [--solver explicit|implicit] [--steps N] [--disp value]
                                  [--backend auto|metal|cpu]
                                  [--nx N] [--ny N] [--order linear|quadratic] [--subdivide L]
                                  [--benchmarks]
                                  [--explicit-substeps N] [--relax-iters N] [--dt value] [--damping value]
                                  [--mass-density value] [--density-penalty value] [--velocity-clamp value]
                                  [--topopt] [--topopt-iters N] [--patch-radius N]
                                  [--min-density value] [--max-density value] [--volume-fraction value]
                                  [--move-limit value] [--reference-stride N]
                                  [--objective-tol value] [--density-change-tol value]
                                  [--objective compliance|mean_von_mises|max_damage|compliance_plus_von_mises]
                                  [--objective-weight value]
                                  [--topopt-export-every N]
                                  [--visualize output-dir] [--deformation-scale value]
                                  [--image-width W] [--image-height H]

            Examples:
              swift run mayor-fem --solver explicit --steps 12 --disp 0.08 --backend metal
              swift run mayor-fem --benchmarks --solver explicit --backend metal
              swift run mayor-fem --topopt --topopt-iters 6 --patch-radius 1 --objective compliance --volume-fraction 0.4
              swift run mayor-fem --visualize out/viz2d --deformation-scale 12
            """
        )
    }
}

private func sumXReactionsOnFixedFace(problem: FEMProblem2D, reactions: [Float]) -> Float {
    var total: Float = 0
    for bc in problem.prescribedDisplacements where bc.component == 0 {
        if abs(bc.value) < 1e-8 {
            total += reactions[bc.dof]
        }
    }
    return total
}

private func writeDensityCSV(_ densities: [Float], outputDirectory: String) throws {
    let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: directoryURL, withIntermediateDirectories: true)
    let fileURL = directoryURL.appendingPathComponent("densities.csv")

    var lines = ["element_id,density"]
    lines.reserveCapacity(densities.count + 1)
    for (index, density) in densities.enumerated() {
        lines.append("\(index),\(density)")
    }
    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

private func writeDensityHistoryCSVs(_ densityHistory: [[Float]], outputDirectory: String) throws {
    let root = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
    let historyDirectory = root.appendingPathComponent("density_history", isDirectory: true)
    try FileManager.default.createDirectory(at: historyDirectory, withIntermediateDirectories: true)

    for (iteration, densities) in densityHistory.enumerated() {
        let fileURL = historyDirectory.appendingPathComponent(
            String(format: "densities_iter_%04d.csv", iteration)
        )
        var lines = ["element_id,density"]
        lines.reserveCapacity(densities.count + 1)
        for (index, density) in densities.enumerated() {
            lines.append("\(index),\(density)")
        }
        try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
    }
}

private func writeTopologyHistoryCSV(_ history: [TopologyIterationResult2D], outputDirectory: String) throws {
    let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: directoryURL, withIntermediateDirectories: true)
    let fileURL = directoryURL.appendingPathComponent("topopt_history.csv")

    var lines = ["iteration,objective,average_density,total_density_change,volume_violation,updated_element_count"]
    lines.reserveCapacity(history.count + 1)
    for item in history {
        lines.append(
            "\(item.iteration),\(item.objective),\(item.averageDensity),\(item.totalDensityChange),\(item.volumeViolation),\(item.updatedElementCount)"
        )
    }
    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

do {
    let options = try CLIOptions.parse(arguments: CommandLine.arguments)

    if options.runBenchmarks {
        let results: [BenchmarkResult2D]
        switch options.solverMode {
        case .implicitNewton:
            results = try LiteratureBenchmarks2D.runAll(backendChoice: options.backend)
        case .explicitDynamics:
            results = try ExplicitLiteratureBenchmarks2D.runAll(backendChoice: options.backend)
        }

        print("2D benchmark suite (\(options.solverMode == .explicitDynamics ? "explicit" : "implicit")):")
        var allPassed = true
        for result in results {
            let label = result.passed ? "PASS" : "FAIL"
            print("  [\(label)] \(result.name)")
            print("         \(result.detail)")
            allPassed = allPassed && result.passed
        }

        if !allPassed {
            exit(2)
        }
        exit(0)
    }

    let problem = ExampleProblems2D.displacementControlledTension(
        nx: options.nx,
        ny: options.ny,
        order: options.order,
        subdivisionLevels: options.subdivisionLevels,
        endDisplacement: options.displacement,
        loadSteps: options.steps
    )

    let result: SolveResult2D
    var topologyResult: TopologyOptimizationResult2D?

    if options.runTopologyOptimization {
        if options.solverMode != .explicitDynamics {
            throw FEMError.invalidBoundaryCondition("Topology optimization currently requires --solver explicit.")
        }
        let objective: PatchObjectiveFunction2D
        switch options.objectiveName {
        case "compliance":
            objective = TopologyObjectives2D.compliance
        case "mean_von_mises":
            objective = TopologyObjectives2D.meanVonMises
        case "max_damage":
            objective = TopologyObjectives2D.maxDamage
        case "compliance_plus_von_mises":
            objective = TopologyObjectives2D.compliancePlusMeanVonMises(weight: options.objectiveWeight)
        default:
            throw FEMError.invalidBoundaryCondition(
                "--objective supports: compliance|mean_von_mises|max_damage|compliance_plus_von_mises"
            )
        }

        let optimizer = try PatchTopologyOptimizer2D(
            problem: problem,
            backendChoice: options.backend,
            controls: options.topologyControls(),
            objective: objective
        )
        let topopt = try optimizer.optimize()
        result = topopt.finalSolve
        topologyResult = topopt

        print("Topology optimization history:")
        for item in topopt.history {
            print(
                String(
                    format: "  iter %2d | objective=%.6e | avg_density=%.4f | total_change=%.4f | vol_violation=%.4e | updates=%d",
                    item.iteration,
                    item.objective,
                    item.averageDensity,
                    item.totalDensityChange,
                    item.volumeViolation,
                    item.updatedElementCount
                )
            )
        }
        print("Topology converged: \(topopt.converged)")
        print("Topology convergence reason: \(topopt.convergenceReason)")
    } else {
        switch options.solverMode {
        case .implicitNewton:
            let solver = try NonlinearFEMSolver2D(problem: problem, backendChoice: options.backend)
            result = try solver.solve()
        case .explicitDynamics:
            let solver = try ExplicitFEMSolver2D(
                problem: problem,
                explicitControls: options.explicitControls(),
                backendChoice: options.backend
            )
            result = try solver.solve()
        }
    }

    print("Backend: \(result.backendName)")
    print("Converged: \(result.converged)")
    print("Step history:")

    for step in result.stepHistory {
        print(
            String(
                format: "  step %2d | load=%.3f | iters=%2d | res=%.3e | eqp_max=%.4f | dmg_max=%.4f",
                step.step,
                step.loadFactor,
                step.iterations,
                step.residualNorm,
                step.maxEquivalentPlasticStrain,
                step.maxDamage
            )
        )
    }

    let supportReaction = sumXReactionsOnFixedFace(problem: problem, reactions: result.reactions)
    let finalDisp = result.displacements.map(\.x).max() ?? 0
    print(String(format: "Final prescribed displacement: %.5f", finalDisp))
    print(String(format: "Support reaction (x, fixed face sum): %.5f", supportReaction))

    if let visualizationDirectory = options.visualizationDirectory {
        let vtkDirectory = URL(fileURLWithPath: visualizationDirectory, isDirectory: true)
            .appendingPathComponent("vtk")
        let pngDirectory = URL(fileURLWithPath: visualizationDirectory, isDirectory: true)
            .appendingPathComponent("png")

        let vtkFiles = try FEMVisualization.writeVTKSeries(
            problem: problem,
            result: result,
            outputDirectory: vtkDirectory.path,
            deformationScale: options.deformationScale
        )

        let pngFiles = try FEMVisualization.writePNGSeries(
            problem: problem,
            result: result,
            outputDirectory: pngDirectory.path,
            deformationScale: options.deformationScale,
            imageWidth: options.imageWidth,
            imageHeight: options.imageHeight
        )

        try writeDensityCSV(result.elementDensities, outputDirectory: visualizationDirectory)
        if let topopt = topologyResult {
            try writeTopologyHistoryCSV(topopt.history, outputDirectory: visualizationDirectory)
            try writeDensityHistoryCSVs(topopt.densityHistory, outputDirectory: visualizationDirectory)

            if options.topologyExportEvery > 0 {
                let iterationsRoot = URL(fileURLWithPath: visualizationDirectory, isDirectory: true)
                    .appendingPathComponent("topopt_iterations", isDirectory: true)
                try FileManager.default.createDirectory(at: iterationsRoot, withIntermediateDirectories: true)

                for (iteration, densities) in topopt.densityHistory.enumerated() {
                    let isFinal = iteration == topopt.densityHistory.count - 1
                    if !isFinal && iteration % options.topologyExportEvery != 0 {
                        continue
                    }

                    let solver = try ExplicitFEMSolver2D(
                        problem: problem,
                        explicitControls: options.explicitControls(),
                        backendChoice: options.backend
                    )
                    let iterationResult = try solver.solve(
                        densities: densities,
                        minimumDensity: options.topologyMinimumDensity
                    )
                    let iterationDir = iterationsRoot.appendingPathComponent(
                        String(format: "iter_%04d", iteration),
                        isDirectory: true
                    )
                    let iterationVtkDir = iterationDir.appendingPathComponent("vtk", isDirectory: true)
                    let iterationPngDir = iterationDir.appendingPathComponent("png", isDirectory: true)

                    _ = try FEMVisualization.writeVTKSeries(
                        problem: problem,
                        result: iterationResult,
                        outputDirectory: iterationVtkDir.path,
                        deformationScale: options.deformationScale
                    )
                    _ = try FEMVisualization.writePNGSeries(
                        problem: problem,
                        result: iterationResult,
                        outputDirectory: iterationPngDir.path,
                        deformationScale: options.deformationScale,
                        imageWidth: options.imageWidth,
                        imageHeight: options.imageHeight
                    )
                    try writeDensityCSV(iterationResult.elementDensities, outputDirectory: iterationDir.path)
                }
            }
        }

        print("Visualization VTK files written: \(vtkFiles.count)")
        print("Visualization PNG files written: \(pngFiles.count)")
        print("Visualization directory: \(visualizationDirectory)")
    }
} catch {
    fputs("Error: \(error)\n", stderr)
    CLIOptions.printUsage()
    exit(1)
}
