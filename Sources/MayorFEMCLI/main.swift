import Foundation
import MayorFEM

private struct CLIOptions {
    var steps: Int = 10
    var displacement: Float = 0.08
    var backend: ComputeBackendChoice = .auto
    var runBenchmarks: Bool = false
    var visualizationDirectory: String?
    var deformationScale: Float = 8.0
    var imageWidth: Int = 1400
    var imageHeight: Int = 900

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
                case "auto":
                    options.backend = .auto
                case "metal":
                    options.backend = .metal
                case "cpu":
                    options.backend = .cpu
                default:
                    throw FEMError.invalidBoundaryCondition("--backend requires one of auto|metal|cpu")
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
                    throw FEMError.invalidBoundaryCondition("--deformation-scale requires a positive floating point value")
                }
                options.deformationScale = value
            case "--image-width":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 256 else {
                    throw FEMError.invalidBoundaryCondition("--image-width requires an integer >= 256")
                }
                options.imageWidth = value
            case "--image-height":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 256 else {
                    throw FEMError.invalidBoundaryCondition("--image-height requires an integer >= 256")
                }
                options.imageHeight = value
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

    static func printUsage() {
        print("""
        mayor-fem

        Usage:
          swift run mayor-fem [--steps N] [--disp value] [--backend auto|metal|cpu] [--benchmarks]
                              [--visualize output-dir] [--deformation-scale value]
                              [--image-width W] [--image-height H]

        Example:
          swift run mayor-fem --steps 12 --disp 0.08 --backend auto
          swift run mayor-fem --benchmarks --backend cpu
          swift run mayor-fem --steps 12 --visualize out/viz --deformation-scale 10
        """)
    }
}

func sumXReactionsOnFixedFace(problem: FEMProblem, reactions: [Float]) -> Float {
    var total: Float = 0
    for bc in problem.prescribedDisplacements where bc.component == 0 {
        if abs(problem.mesh.nodes[bc.node].x) < 1e-6 {
            total += reactions[bc.dof]
        }
    }
    return total
}

do {
    let options = try CLIOptions.parse(arguments: CommandLine.arguments)

    if options.runBenchmarks {
        let results = try LiteratureBenchmarks.runAll(backendChoice: options.backend)
        print("Benchmark suite:")

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

    let problem = ExampleProblems.displacementControlledTension(
        endDisplacement: options.displacement,
        loadSteps: options.steps
    )

    let solver = try NonlinearFEMSolver(problem: problem, backendChoice: options.backend)
    let result = try solver.solve()

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

        print("Visualization VTK files written: \(vtkFiles.count)")
        print("Visualization PNG files written: \(pngFiles.count)")
        print("Visualization directory: \(visualizationDirectory)")
    }
} catch {
    fputs("Error: \(error)\n", stderr)
    CLIOptions.printUsage()
    exit(1)
}
