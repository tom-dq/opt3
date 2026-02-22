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
    var integrationScheme: IntegrationScheme2D = .full
    var subdivisionLevels: Int = 0

    var explicitSubsteps: Int = 48
    var explicitRelaxIterations: Int = 6
    var explicitTimeStep: Float = 4e-4
    var explicitDamping: Float = 0.08
    var explicitMassDensity: Float = 1_250.0
    var explicitDensityPenalty: Float = 3.0
    var explicitVelocityClamp: Float = 25.0
    var stabilizationResidualBlend: Float = 0.0
    var stabilizationVelocitySmoothing: Float = 0.0
    var stabilizationDisplacementSmoothing: Float = 0.0
    var stabilizationSmoothingPasses: Int = 1

    var runTopologyOptimization: Bool = false
    var runStabilizationBenchmark: Bool = false
    var runTopologyLiteratureBenchmark: Bool = false
    var topologyIterations: Int = 8
    var topologyPatchRadius: Int = 1
    var topologyMinimumDensity: Float = 0.001
    var topologyMaximumDensity: Float = 1.0
    var topologyTargetVolumeFractionStart: Float?
    var topologyTargetVolumeFractionEnd: Float?
    var topologyMoveLimit: Float = 0.12
    var topologyReferenceStride: Int = 1
    var topologyMaxReferenceEvaluations: Int?
    var topologyObjectiveTolerance: Float = 1e-5
    var topologyDensityChangeTolerance: Float = 1e-4
    var topologyExportEvery: Int = 0
    var objectiveName: String = "compliance"
    var objectiveWeight: Float = 1e-3
    var convoyWorkers: Int = max(1, ProcessInfo.processInfo.activeProcessorCount / 2)
    var literatureCaseFilter: Set<String> = []
    var literatureLoads: [Float] = [5, 10, 15, 20]
    var literatureLBracketNX: Int = 24
    var literatureLBracketNY: Int = 24
    var literatureUBracketNX: Int = 28
    var literatureUBracketNY: Int = 20
    var literatureCantileverNX: Int = 36
    var literatureCantileverNY: Int = 12
    var literatureStressWeight: Float = 0.05
    var literaturePlasticWeight: Float = 0
    var literatureDamageWeight: Float = 0
    var literatureDensityVarianceWeight: Float = 12

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
            case "--integration":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--integration requires full|reduced")
                }
                switch arguments[index] {
                case "full": options.integrationScheme = .full
                case "reduced": options.integrationScheme = .reduced
                default:
                    throw FEMError.invalidBoundaryCondition("--integration requires full|reduced")
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
            case "--residual-blend":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0, value < 1 else {
                    throw FEMError.invalidBoundaryCondition("--residual-blend requires 0 <= value < 1")
                }
                options.stabilizationResidualBlend = value
            case "--velocity-smoothing":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0, value < 1 else {
                    throw FEMError.invalidBoundaryCondition("--velocity-smoothing requires 0 <= value < 1")
                }
                options.stabilizationVelocitySmoothing = value
            case "--displacement-smoothing":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0, value < 1 else {
                    throw FEMError.invalidBoundaryCondition("--displacement-smoothing requires 0 <= value < 1")
                }
                options.stabilizationDisplacementSmoothing = value
            case "--smoothing-passes":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--smoothing-passes requires integer > 0")
                }
                options.stabilizationSmoothingPasses = value
            case "--stabilization-bench":
                options.runStabilizationBenchmark = true
            case "--topopt-literature-bench":
                options.runTopologyLiteratureBenchmark = true
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
                options.topologyTargetVolumeFractionStart = value
                options.topologyTargetVolumeFractionEnd = value
            case "--volume-fraction-start":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0, value <= 1 else {
                    throw FEMError.invalidBoundaryCondition("--volume-fraction-start requires 0 < value <= 1")
                }
                options.topologyTargetVolumeFractionStart = value
            case "--volume-fraction-end":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value > 0, value <= 1 else {
                    throw FEMError.invalidBoundaryCondition("--volume-fraction-end requires 0 < value <= 1")
                }
                options.topologyTargetVolumeFractionEnd = value
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
            case "--max-reference-evals":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--max-reference-evals requires integer > 0")
                }
                options.topologyMaxReferenceEvaluations = value
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
            case "--convoy-workers":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value > 0 else {
                    throw FEMError.invalidBoundaryCondition("--convoy-workers requires integer > 0")
                }
                options.convoyWorkers = value
            case "--literature-cases":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--literature-cases requires comma-separated ids")
                }
                let ids = arguments[index]
                    .split(separator: ",")
                    .map { String($0).trimmingCharacters(in: .whitespacesAndNewlines) }
                    .filter { !$0.isEmpty }
                options.literatureCaseFilter = Set(ids)
            case "--literature-loads":
                index += 1
                guard index < arguments.count else {
                    throw FEMError.invalidBoundaryCondition("--literature-loads requires comma-separated positive values")
                }
                let parsedLoads = try arguments[index]
                    .split(separator: ",")
                    .map { raw -> Float in
                        guard let value = Float(raw.trimmingCharacters(in: .whitespacesAndNewlines)), value > 0 else {
                            throw FEMError.invalidBoundaryCondition("--literature-loads requires comma-separated positive values")
                        }
                        return value
                    }
                if !parsedLoads.isEmpty {
                    options.literatureLoads = parsedLoads
                }
            case "--literature-lbracket-nx":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 8 else {
                    throw FEMError.invalidBoundaryCondition("--literature-lbracket-nx requires integer >= 8")
                }
                options.literatureLBracketNX = value
            case "--literature-lbracket-ny":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 8 else {
                    throw FEMError.invalidBoundaryCondition("--literature-lbracket-ny requires integer >= 8")
                }
                options.literatureLBracketNY = value
            case "--literature-ubracket-nx":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 8 else {
                    throw FEMError.invalidBoundaryCondition("--literature-ubracket-nx requires integer >= 8")
                }
                options.literatureUBracketNX = value
            case "--literature-ubracket-ny":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 8 else {
                    throw FEMError.invalidBoundaryCondition("--literature-ubracket-ny requires integer >= 8")
                }
                options.literatureUBracketNY = value
            case "--literature-cantilever-nx":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 8 else {
                    throw FEMError.invalidBoundaryCondition("--literature-cantilever-nx requires integer >= 8")
                }
                options.literatureCantileverNX = value
            case "--literature-cantilever-ny":
                index += 1
                guard index < arguments.count, let value = Int(arguments[index]), value >= 4 else {
                    throw FEMError.invalidBoundaryCondition("--literature-cantilever-ny requires integer >= 4")
                }
                options.literatureCantileverNY = value
            case "--literature-stress-weight":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--literature-stress-weight requires value >= 0")
                }
                options.literatureStressWeight = value
            case "--literature-plastic-weight":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--literature-plastic-weight requires value >= 0")
                }
                options.literaturePlasticWeight = value
            case "--literature-damage-weight":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--literature-damage-weight requires value >= 0")
                }
                options.literatureDamageWeight = value
            case "--literature-density-variance-weight":
                index += 1
                guard index < arguments.count, let value = Float(arguments[index]), value >= 0 else {
                    throw FEMError.invalidBoundaryCondition("--literature-density-variance-weight requires value >= 0")
                }
                options.literatureDensityVarianceWeight = value
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

    func stabilizationControls() -> ExplicitStabilizationControls2D {
        ExplicitStabilizationControls2D(
            residualBlend: stabilizationResidualBlend,
            velocitySmoothing: stabilizationVelocitySmoothing,
            displacementSmoothing: stabilizationDisplacementSmoothing,
            smoothingPasses: stabilizationSmoothingPasses
        )
    }

    func explicitControls() -> ExplicitSolverControls2D {
        ExplicitSolverControls2D(
            stabilization: stabilizationControls(),
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
            targetVolumeFractionStart: topologyTargetVolumeFractionStart,
            targetVolumeFractionEnd: topologyTargetVolumeFractionEnd,
            moveLimit: topologyMoveLimit,
            referenceElementStride: topologyReferenceStride,
            maxReferenceEvaluationsPerIteration: topologyMaxReferenceEvaluations,
            objectiveTolerance: topologyObjectiveTolerance,
            densityChangeTolerance: topologyDensityChangeTolerance,
            explicitControls: explicitControls()
        )
    }

    func topologyLiteratureControls() -> TopologyLiteratureBenchmarkControls2D {
        let literatureIterations = topologyIterations == 8 ? 4 : topologyIterations
        let literatureStride = topologyReferenceStride == 1 ? 12 : topologyReferenceStride
        return TopologyLiteratureBenchmarkControls2D(
            order: order,
            integrationScheme: integrationScheme,
            subdivisionLevels: subdivisionLevels,
            iterations: literatureIterations,
            patchRadius: topologyPatchRadius,
            minimumDensity: topologyMinimumDensity,
            maximumDensity: topologyMaximumDensity,
            moveLimit: topologyMoveLimit,
            referenceElementStride: literatureStride,
            maxReferenceEvaluationsPerIteration: topologyMaxReferenceEvaluations ?? 24,
            objectiveTolerance: topologyObjectiveTolerance,
            densityChangeTolerance: topologyDensityChangeTolerance,
            targetVolumeFractionStart: topologyTargetVolumeFractionStart,
            targetVolumeFractionEnd: topologyTargetVolumeFractionEnd,
            explicitControls: explicitControls(),
            convoyWorkers: convoyWorkers,
            caseIDs: literatureCaseFilter.isEmpty ? nil : literatureCaseFilter,
            liangLoads: literatureLoads,
            lBracketNX: literatureLBracketNX,
            lBracketNY: literatureLBracketNY,
            uBracketNX: literatureUBracketNX,
            uBracketNY: literatureUBracketNY,
            cantileverNX: literatureCantileverNX,
            cantileverNY: literatureCantileverNY,
            loadSteps: steps,
            stressWeight: literatureStressWeight,
            plasticWeight: literaturePlasticWeight,
            damageWeight: literatureDamageWeight,
            densityVarianceWeight: literatureDensityVarianceWeight
        )
    }

    static func printUsage() {
        print(
            """
            mayor-fem (2D)

            Usage:
              swift run mayor-fem [--solver explicit|implicit] [--steps N] [--disp value]
                                  [--backend auto|metal|cpu]
                                  [--nx N] [--ny N] [--order linear|quadratic]
                                  [--integration full|reduced] [--subdivide L]
                                  [--benchmarks]
                                  [--explicit-substeps N] [--relax-iters N] [--dt value] [--damping value]
                                  [--mass-density value] [--density-penalty value] [--velocity-clamp value]
                                  [--residual-blend value] [--velocity-smoothing value]
                                  [--displacement-smoothing value] [--smoothing-passes N]
                                  [--stabilization-bench]
                                  [--topopt-literature-bench] [--convoy-workers N]
                                  [--literature-cases id1,id2] [--literature-loads l1,l2,...]
                                  [--literature-lbracket-nx N] [--literature-lbracket-ny N]
                                  [--literature-ubracket-nx N] [--literature-ubracket-ny N]
                                  [--literature-cantilever-nx N] [--literature-cantilever-ny N]
                                  [--literature-stress-weight value]
                                  [--literature-plastic-weight value]
                                  [--literature-damage-weight value]
                                  [--literature-density-variance-weight value]
                                  [--topopt] [--topopt-iters N] [--patch-radius N]
                                  [--min-density value] [--max-density value]
                                  [--volume-fraction value]
                                  [--volume-fraction-start value] [--volume-fraction-end value]
                                  [--move-limit value] [--reference-stride N]
                                  [--max-reference-evals N]
                                  [--objective-tol value] [--density-change-tol value]
                                  [--objective compliance|mean_von_mises|max_damage|compliance_plus_von_mises|compliance_plus_plasticity]
                                  [--objective-weight value]
                                  [--topopt-export-every N]
                                  [--visualize output-dir] [--deformation-scale value]
                                  [--image-width W] [--image-height H]

            Examples:
              swift run mayor-fem --solver explicit --steps 12 --disp 0.08 --backend metal
              swift run mayor-fem --benchmarks --solver explicit --backend metal
              swift run mayor-fem --stabilization-bench --backend cpu --visualize out/stabilization_sweep
              swift run mayor-fem --topopt --topopt-iters 6 --patch-radius 1 --objective compliance --volume-fraction 0.4
              swift run mayor-fem --topopt --topopt-iters 12 --volume-fraction-start 0.8 --volume-fraction-end 0.35
              swift run mayor-fem --topopt-literature-bench --solver explicit --backend metal --convoy-workers 3 --visualize out/topopt_literature
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

    var lines = ["iteration,objective,average_density,target_volume_fraction,total_density_change,volume_violation,updated_element_count"]
    lines.reserveCapacity(history.count + 1)
    for item in history {
        lines.append(
            "\(item.iteration),\(item.objective),\(item.averageDensity),\(item.targetVolumeFraction),\(item.totalDensityChange),\(item.volumeViolation),\(item.updatedElementCount)"
        )
    }
    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

private func writeStabilizationBenchmarkCSV(_ sweep: StabilizationBenchmarkSweep2D, outputDirectory: String) throws {
    let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: directoryURL, withIntermediateDirectories: true)
    let fileURL = directoryURL.appendingPathComponent("stabilization_benchmark.csv")

    var lines = [
        "case,profile,roughness_index,refinement_mismatch,final_residual_norm,max_displacement,runtime_seconds,roughness_reduction_percent,mismatch_reduction_percent,combined_effectiveness_percent"
    ]
    lines.reserveCapacity(sweep.rows.count + 1)
    for row in sweep.rows {
        lines.append(
            "\(row.caseName),\(row.profileName),\(row.roughnessIndex),\(row.refinementMismatch),\(row.finalResidualNorm),\(row.maxDisplacementMagnitude),\(row.runtimeSeconds),\(row.roughnessReductionPercent),\(row.mismatchReductionPercent),\(row.combinedEffectivenessPercent)"
        )
    }
    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

private func printStabilizationBenchmarkTable(_ sweep: StabilizationBenchmarkSweep2D) {
    print("| Case | Profile | Roughness | Refinement Mismatch | Residual | Max | Runtime (s) | Roughness Δ% | Mismatch Δ% | Combined Δ% |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in sweep.rows {
        print(
            String(
                format: "| %@ | %@ | %.6f | %.6f | %.4e | %.5f | %.2f | %.2f | %.2f | %.2f |",
                row.caseName,
                row.profileName,
                row.roughnessIndex,
                row.refinementMismatch,
                row.finalResidualNorm,
                row.maxDisplacementMagnitude,
                row.runtimeSeconds,
                row.roughnessReductionPercent,
                row.mismatchReductionPercent,
                row.combinedEffectivenessPercent
            )
        )
    }
}

private func writeTopologyLiteratureRunsCSV(
    _ runs: [TopologyLiteratureRun2D],
    outputDirectory: String
) throws {
    let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: directoryURL, withIntermediateDirectories: true)
    let fileURL = directoryURL.appendingPathComponent("topopt_literature_runs.csv")

    var lines = [
        "case_id,case_name,reference,profile,load_magnitude,element_count,runtime_seconds,objective,compliance,average_density,final_residual,max_von_mises,max_eq_plastic,max_damage,yield_onset_load_factor,right_edge_density,converged,iterations_completed"
    ]
    lines.reserveCapacity(runs.count + 1)

    for run in runs {
        let loadMagnitude = run.loadMagnitude.map { String($0) } ?? ""
        lines.append(
            "\"\(run.caseID)\",\"\(run.caseName)\",\"\(run.reference)\",\"\(run.profileName)\",\(loadMagnitude),\(run.elementCount),\(run.runtimeSeconds),\(run.objective),\(run.compliance),\(run.averageDensity),\(run.finalResidualNorm),\(run.maxVonMises),\(run.maxEquivalentPlasticStrain),\(run.maxDamage),\(run.yieldOnsetLoadFactor),\(run.rightEdgeDensity),\(run.converged),\(run.iterationsCompleted)"
        )
    }

    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

private func writeTopologyLiteratureChecksCSV(
    _ checks: [TopologyLiteratureTrendCheck2D],
    outputDirectory: String
) throws {
    let directoryURL = URL(fileURLWithPath: outputDirectory, isDirectory: true)
    try FileManager.default.createDirectory(at: directoryURL, withIntermediateDirectories: true)
    let fileURL = directoryURL.appendingPathComponent("topopt_literature_checks.csv")

    var lines = ["case_name,reference,metric,observed,target,passed"]
    lines.reserveCapacity(checks.count + 1)
    for check in checks {
        lines.append(
            "\"\(check.caseName)\",\"\(check.reference)\",\"\(check.metric)\",\"\(check.observed)\",\"\(check.target)\",\(check.passed)"
        )
    }

    try lines.joined(separator: "\n").appending("\n").write(to: fileURL, atomically: true, encoding: .utf8)
}

private func printTopologyLiteratureRunsTable(_ runs: [TopologyLiteratureRun2D]) {
    print("| Case | Profile | Load | Elements | Compliance | Avg Density | Max EqP | Yield LF | Right-Edge ρ | Runtime (s) |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")

    for run in runs {
        let load = run.loadMagnitude.map { String(format: "%.1f", $0) } ?? "-"
        print(
            String(
                format: "| %@ | %@ | %@ | %d | %.6e | %.4f | %.6f | %.3f | %.4f | %.2f |",
                run.caseName,
                run.profileName,
                load,
                run.elementCount,
                run.compliance,
                run.averageDensity,
                run.maxEquivalentPlasticStrain,
                run.yieldOnsetLoadFactor,
                run.rightEdgeDensity,
                run.runtimeSeconds
            )
        )
    }
}

private func printTopologyLiteratureChecksTable(_ checks: [TopologyLiteratureTrendCheck2D]) {
    print("| Case | Metric | Observed | Target | Status |")
    print("|---|---|---:|---:|---|")
    for check in checks {
        print(
            "| \(check.caseName) | \(check.metric) | \(check.observed) | \(check.target) | \(check.passed ? "PASS" : "FAIL") |"
        )
    }
}

private func exportTopologyLiteratureVisuals(
    runs: [TopologyLiteratureRun2D],
    outputDirectory: String,
    deformationScale: Float,
    imageWidth: Int,
    imageHeight: Int
) throws {
    let root = URL(fileURLWithPath: outputDirectory, isDirectory: true)
        .appendingPathComponent("literature_runs", isDirectory: true)
    try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)

    for run in runs {
        let loadLabel = run.loadMagnitude.map { String(format: "_load_%04.1f", $0).replacingOccurrences(of: ".", with: "p") } ?? ""
        let runID = "\(run.caseID)__\(run.profileName)\(loadLabel)"
        let runDirectory = root.appendingPathComponent(runID, isDirectory: true)
        let vtkDirectory = runDirectory.appendingPathComponent("vtk", isDirectory: true)
        let pngDirectory = runDirectory.appendingPathComponent("png", isDirectory: true)

        _ = try FEMVisualization.writeVTKSeries(
            problem: run.problem,
            result: run.topologyResult.finalSolve,
            outputDirectory: vtkDirectory.path,
            deformationScale: deformationScale
        )
        _ = try FEMVisualization.writePNGSeries(
            problem: run.problem,
            result: run.topologyResult.finalSolve,
            outputDirectory: pngDirectory.path,
            deformationScale: deformationScale,
            imageWidth: imageWidth,
            imageHeight: imageHeight
        )

        try writeDensityCSV(run.topologyResult.densities, outputDirectory: runDirectory.path)
        try writeTopologyHistoryCSV(run.topologyResult.history, outputDirectory: runDirectory.path)
        try writeDensityHistoryCSVs(run.topologyResult.densityHistory, outputDirectory: runDirectory.path)
    }
}

do {
    let options = try CLIOptions.parse(arguments: CommandLine.arguments)

    if options.runStabilizationBenchmark {
        guard options.solverMode == .explicitDynamics else {
            throw FEMError.invalidBoundaryCondition("--stabilization-bench requires --solver explicit.")
        }

        let sweep = try StabilizationBenchmarks2D.runSweep(backendChoice: options.backend)
        print("Stabilization benchmark sweep (lower roughness/mismatch is better):")
        printStabilizationBenchmarkTable(sweep)
        if let visualizationDirectory = options.visualizationDirectory {
            try writeStabilizationBenchmarkCSV(sweep, outputDirectory: visualizationDirectory)
            print("Stabilization benchmark CSV: \(visualizationDirectory)/stabilization_benchmark.csv")
        }
        exit(0)
    }

    if options.runTopologyLiteratureBenchmark {
        guard options.solverMode == .explicitDynamics else {
            throw FEMError.invalidBoundaryCondition("--topopt-literature-bench requires --solver explicit.")
        }

        let sweep = try TopologyLiteratureBenchmarks2D.runConvoy(
            backendChoice: options.backend,
            controls: options.topologyLiteratureControls()
        )

        print("Topology literature convoy runs:")
        printTopologyLiteratureRunsTable(sweep.runs)
        print("")
        print("Topology literature trend checks:")
        printTopologyLiteratureChecksTable(sweep.checks)

        if let visualizationDirectory = options.visualizationDirectory {
            try writeTopologyLiteratureRunsCSV(sweep.runs, outputDirectory: visualizationDirectory)
            try writeTopologyLiteratureChecksCSV(sweep.checks, outputDirectory: visualizationDirectory)
            try exportTopologyLiteratureVisuals(
                runs: sweep.runs,
                outputDirectory: visualizationDirectory,
                deformationScale: options.deformationScale,
                imageWidth: options.imageWidth,
                imageHeight: options.imageHeight
            )
            print("Topology literature CSV: \(visualizationDirectory)/topopt_literature_runs.csv")
            print("Topology literature checks CSV: \(visualizationDirectory)/topopt_literature_checks.csv")
            print("Topology literature visualizations: \(visualizationDirectory)/literature_runs/")
        }

        if !sweep.allPassed {
            exit(2)
        }
        exit(0)
    }

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
        integrationScheme: options.integrationScheme,
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
        case "compliance_plus_plasticity":
            objective = TopologyObjectives2D.compliancePlusPlasticity(
                stressWeight: options.objectiveWeight,
                plasticWeight: 250 * options.objectiveWeight,
                damageWeight: 80 * options.objectiveWeight
            )
        default:
            throw FEMError.invalidBoundaryCondition(
                "--objective supports: compliance|mean_von_mises|max_damage|compliance_plus_von_mises|compliance_plus_plasticity"
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
                    format: "  iter %2d | objective=%.6e | avg_density=%.4f | target=%.4f | total_change=%.4f | vol_violation=%.4e | updates=%d",
                    item.iteration,
                    item.objective,
                    item.averageDensity,
                    item.targetVolumeFraction,
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
