import Foundation

public enum TopologyObjectives2D {
    public static func compliance(_ input: PatchObjectiveInput2D) -> Float {
        input.patchElements.reduce(Float(0)) { partial, elementIndex in
            guard elementIndex >= 0, elementIndex < input.solveResult.elementEnergies.count else {
                return partial
            }
            return partial + input.solveResult.elementEnergies[elementIndex]
        }
    }
}

public final class PatchTopologyOptimizer2D {
    private let problem: FEMProblem2D
    private let backendChoice: ComputeBackendChoice
    private let controls: TopologyOptimizationControls2D
    private let objective: PatchObjectiveFunction2D
    private let elementAdjacency: [[Int]]
    private let elementCount: Int

    public init(
        problem: FEMProblem2D,
        backendChoice: ComputeBackendChoice = .auto,
        controls: TopologyOptimizationControls2D = TopologyOptimizationControls2D(),
        objective: @escaping PatchObjectiveFunction2D = TopologyObjectives2D.compliance
    ) throws {
        self.problem = problem
        self.backendChoice = backendChoice
        self.controls = controls
        self.objective = objective

        let preparedMesh = try problem.mesh.prepare()
        self.elementCount = preparedMesh.elements.count
        self.elementAdjacency = PatchTopologyOptimizer2D.buildElementAdjacency(preparedMesh: preparedMesh)
    }

    public func optimize(initialDensities: [Float]? = nil) throws -> TopologyOptimizationResult2D {
        var densities = initialDensities ?? Array(repeating: Float(1), count: elementCount)
        if densities.count != elementCount {
            throw FEMError.invalidMesh("Topology optimization density count does not match element count.")
        }
        for index in densities.indices {
            densities[index] = max(controls.minimumDensity, min(controls.maximumDensity, densities[index]))
        }

        var history: [TopologyIterationResult2D] = []
        var finalSolve = try solveForDensities(densities)

        for iteration in 1...controls.iterations {
            var updated = densities
            var totalDensityChange: Float = 0

            for referenceElement in stride(from: 0, to: elementCount, by: controls.referenceElementStride) {
                let patch = patchElements(around: referenceElement, radius: controls.patchRadius)
                if patch.isEmpty {
                    continue
                }

                var onDensities = densities
                onDensities[referenceElement] = controls.maximumDensity
                let onSolve = try solveForDensities(onDensities)
                let onObjective = objective(
                    PatchObjectiveInput2D(
                        problem: problem,
                        solveResult: onSolve,
                        patchElements: patch,
                        densities: onDensities,
                        referenceElement: referenceElement
                    )
                )

                var offDensities = densities
                offDensities[referenceElement] = controls.minimumDensity
                let offSolve = try solveForDensities(offDensities)
                let offObjective = objective(
                    PatchObjectiveInput2D(
                        problem: problem,
                        solveResult: offSolve,
                        patchElements: patch,
                        densities: offDensities,
                        referenceElement: referenceElement
                    )
                )

                let current = densities[referenceElement]
                if onObjective < offObjective {
                    let next = min(controls.maximumDensity, current + controls.moveLimit)
                    updated[referenceElement] = next
                    totalDensityChange += abs(next - current)
                } else if offObjective < onObjective {
                    let next = max(controls.minimumDensity, current - controls.moveLimit)
                    updated[referenceElement] = next
                    totalDensityChange += abs(next - current)
                }
            }

            densities = updated
            finalSolve = try solveForDensities(densities)
            let fullPatch = Array(0..<elementCount)
            let objectiveValue = objective(
                PatchObjectiveInput2D(
                    problem: problem,
                    solveResult: finalSolve,
                    patchElements: fullPatch,
                    densities: densities,
                    referenceElement: -1
                )
            )
            let averageDensity = densities.reduce(0, +) / Float(max(1, densities.count))

            history.append(
                TopologyIterationResult2D(
                    iteration: iteration,
                    objective: objectiveValue,
                    averageDensity: averageDensity,
                    totalDensityChange: totalDensityChange
                )
            )
        }

        return TopologyOptimizationResult2D(densities: densities, history: history, finalSolve: finalSolve)
    }

    private func solveForDensities(_ densities: [Float]) throws -> SolveResult2D {
        let solver = try ExplicitFEMSolver2D(
            problem: problem,
            explicitControls: controls.explicitControls,
            backendChoice: backendChoice
        )
        return try solver.solve(densities: densities, minimumDensity: controls.minimumDensity)
    }

    private func patchElements(around referenceElement: Int, radius: Int) -> [Int] {
        if referenceElement < 0 || referenceElement >= elementCount {
            return []
        }
        if radius <= 0 {
            return [referenceElement]
        }

        var visited = Set<Int>([referenceElement])
        var frontier = [referenceElement]

        for _ in 0..<radius {
            var nextFrontier: [Int] = []
            for element in frontier {
                for neighbor in elementAdjacency[element] where !visited.contains(neighbor) {
                    visited.insert(neighbor)
                    nextFrontier.append(neighbor)
                }
            }
            frontier = nextFrontier
            if frontier.isEmpty {
                break
            }
        }

        return visited.sorted()
    }

    private static func buildElementAdjacency(preparedMesh: PreparedMesh2D) -> [[Int]] {
        var nodeToElements = Array(repeating: [Int](), count: preparedMesh.nodes.count)
        for (elementIndex, element) in preparedMesh.elements.enumerated() {
            nodeToElements[Int(element.nodeIDs[0])].append(elementIndex)
            nodeToElements[Int(element.nodeIDs[1])].append(elementIndex)
            nodeToElements[Int(element.nodeIDs[2])].append(elementIndex)
        }

        var adjacency = Array(repeating: Set<Int>(), count: preparedMesh.elements.count)
        for elements in nodeToElements where !elements.isEmpty {
            for element in elements {
                for neighbor in elements where neighbor != element {
                    adjacency[element].insert(neighbor)
                }
            }
        }

        return adjacency.map { $0.sorted() }
    }
}
