import Foundation
import simd

public enum ElementOrder2D {
    case linear
    case quadratic
}

public enum IntegrationScheme2D {
    case full
    case reduced
}

public enum Element2D {
    case quad4(SIMD4<UInt32>)
}

public struct Mesh2D {
    public var nodes: [SIMD2<Float>]
    public var elements: [Element2D]

    public init(nodes: [SIMD2<Float>], elements: [Element2D]) {
        self.nodes = nodes
        self.elements = elements
    }

    public static func rectangularPlate(
        length: Float = 1.0,
        height: Float = 0.25,
        nx: Int = 8,
        ny: Int = 3,
        order: ElementOrder2D = .linear
    ) -> Mesh2D {
        precondition(nx > 0 && ny > 0)

        // Quadratic mode uses an enriched quad mesh with doubled base resolution.
        let effectiveNX = order == .quadratic ? 2 * nx : nx
        let effectiveNY = order == .quadratic ? 2 * ny : ny

        let dx = length / Float(effectiveNX)
        let dy = height / Float(effectiveNY)

        var nodes: [SIMD2<Float>] = []
        nodes.reserveCapacity((effectiveNX + 1) * (effectiveNY + 1))

        for j in 0...effectiveNY {
            for i in 0...effectiveNX {
                nodes.append(SIMD2<Float>(Float(i) * dx, Float(j) * dy))
            }
        }

        func node(_ i: Int, _ j: Int) -> UInt32 {
            UInt32(j * (effectiveNX + 1) + i)
        }

        var elements: [Element2D] = []
        elements.reserveCapacity(effectiveNX * effectiveNY)

        for j in 0..<effectiveNY {
            for i in 0..<effectiveNX {
                let n00 = node(i, j)
                let n10 = node(i + 1, j)
                let n11 = node(i + 1, j + 1)
                let n01 = node(i, j + 1)
                elements.append(.quad4(SIMD4<UInt32>(n00, n10, n11, n01)))
            }
        }

        return Mesh2D(nodes: nodes, elements: elements)
    }

    public func subdivided(levels: Int) -> Mesh2D {
        if levels <= 0 {
            return self
        }

        var workingNodes = nodes
        var quads = linearizedQuads(elements: elements)

        for _ in 0..<levels {
            var edgeMap: [EdgeKey2D: UInt32] = [:]
            var refined: [SIMD4<UInt32>] = []
            refined.reserveCapacity(quads.count * 4)

            func midpoint(_ a: UInt32, _ b: UInt32) -> UInt32 {
                let key = EdgeKey2D(a, b)
                if let id = edgeMap[key] {
                    return id
                }
                let pa = workingNodes[Int(a)]
                let pb = workingNodes[Int(b)]
                let point = 0.5 * (pa + pb)
                let id = UInt32(workingNodes.count)
                workingNodes.append(point)
                edgeMap[key] = id
                return id
            }

            for quad in quads {
                let n0 = quad[0]
                let n1 = quad[1]
                let n2 = quad[2]
                let n3 = quad[3]

                let m01 = midpoint(n0, n1)
                let m12 = midpoint(n1, n2)
                let m23 = midpoint(n2, n3)
                let m30 = midpoint(n3, n0)

                let center = UInt32(workingNodes.count)
                let p0 = workingNodes[Int(n0)]
                let p1 = workingNodes[Int(n1)]
                let p2 = workingNodes[Int(n2)]
                let p3 = workingNodes[Int(n3)]
                workingNodes.append(0.25 * (p0 + p1 + p2 + p3))

                refined.append(SIMD4<UInt32>(n0, m01, center, m30))
                refined.append(SIMD4<UInt32>(m01, n1, m12, center))
                refined.append(SIMD4<UInt32>(center, m12, n2, m23))
                refined.append(SIMD4<UInt32>(m30, center, m23, n3))
            }

            quads = refined
        }

        return Mesh2D(nodes: workingNodes, elements: quads.map(Element2D.quad4))
    }

    func linearizedQuads(elements: [Element2D]) -> [SIMD4<UInt32>] {
        var quads: [SIMD4<UInt32>] = []
        quads.reserveCapacity(elements.count)

        for element in elements {
            switch element {
            case .quad4(let quad):
                quads.append(quad)
            }
        }

        return quads
    }
}

private struct EdgeKey2D: Hashable {
    var a: UInt32
    var b: UInt32

    init(_ x: UInt32, _ y: UInt32) {
        if x < y {
            self.a = x
            self.b = y
        } else {
            self.a = y
            self.b = x
        }
    }
}

public struct PrescribedDisplacement2D {
    public var node: Int
    public var component: Int
    public var value: Float

    public init(node: Int, component: Int, value: Float) {
        self.node = node
        self.component = component
        self.value = value
    }

    public var dof: Int {
        node * 2 + component
    }
}

public struct FEMProblem2D {
    public var mesh: Mesh2D
    public var material: MaterialParameters
    public var thickness: Float
    public var integrationScheme: IntegrationScheme2D
    public var prescribedDisplacements: [PrescribedDisplacement2D]
    public var controls: SolverControls

    public init(
        mesh: Mesh2D,
        material: MaterialParameters,
        thickness: Float = 1.0,
        integrationScheme: IntegrationScheme2D = .full,
        prescribedDisplacements: [PrescribedDisplacement2D],
        controls: SolverControls = SolverControls()
    ) {
        self.mesh = mesh
        self.material = material
        self.thickness = thickness
        self.integrationScheme = integrationScheme
        self.prescribedDisplacements = prescribedDisplacements
        self.controls = controls
    }
}

public struct StepResult2D {
    public var step: Int
    public var loadFactor: Float
    public var iterations: Int
    public var residualNorm: Float
    public var maxEquivalentPlasticStrain: Float
    public var maxDamage: Float
}

public struct StepSnapshot2D {
    public var step: Int
    public var loadFactor: Float
    public var displacements: [SIMD2<Float>]
    public var elementStates: [ElementState]
    public var elementVonMises: [Float]
    public var elementEnergies: [Float]
    public var elementDensities: [Float]
}

public struct SolveResult2D {
    public var displacements: [SIMD2<Float>]
    public var reactions: [Float]
    public var elementStates: [ElementState]
    public var elementVonMises: [Float]
    public var elementEnergies: [Float]
    public var elementDensities: [Float]
    public var stepHistory: [StepResult2D]
    public var stepSnapshots: [StepSnapshot2D]
    public var backendName: String
    public var converged: Bool
}

public struct ExplicitSolverControls2D {
    public var substepsPerLoadStep: Int
    public var relaxationIterationsPerSubstep: Int
    public var timeStep: Float
    public var damping: Float
    public var massDensity: Float
    public var densityPenalty: Float
    public var velocityClamp: Float

    public init(
        substepsPerLoadStep: Int = 48,
        relaxationIterationsPerSubstep: Int = 6,
        timeStep: Float = 4e-4,
        damping: Float = 0.08,
        massDensity: Float = 1_250.0,
        densityPenalty: Float = 3.0,
        velocityClamp: Float = 25.0
    ) {
        self.substepsPerLoadStep = max(1, substepsPerLoadStep)
        self.relaxationIterationsPerSubstep = max(1, relaxationIterationsPerSubstep)
        self.timeStep = max(1e-8, timeStep)
        self.damping = max(0, damping)
        self.massDensity = max(1e-6, massDensity)
        self.densityPenalty = max(1.0, densityPenalty)
        self.velocityClamp = max(1e-4, velocityClamp)
    }
}

public struct TopologyOptimizationControls2D {
    public var iterations: Int
    public var patchRadius: Int
    public var minimumDensity: Float
    public var maximumDensity: Float
    public var targetVolumeFractionStart: Float?
    public var targetVolumeFractionEnd: Float?
    public var moveLimit: Float
    public var referenceElementStride: Int
    public var objectiveTolerance: Float
    public var densityChangeTolerance: Float
    public var explicitControls: ExplicitSolverControls2D

    public init(
        iterations: Int = 8,
        patchRadius: Int = 1,
        minimumDensity: Float = 0.001,
        maximumDensity: Float = 1.0,
        targetVolumeFractionStart: Float? = nil,
        targetVolumeFractionEnd: Float? = nil,
        moveLimit: Float = 0.12,
        referenceElementStride: Int = 1,
        objectiveTolerance: Float = 1e-5,
        densityChangeTolerance: Float = 1e-4,
        explicitControls: ExplicitSolverControls2D = ExplicitSolverControls2D()
    ) {
        self.iterations = max(1, iterations)
        self.patchRadius = max(0, patchRadius)
        self.minimumDensity = max(0.001, min(1.0, minimumDensity))
        self.maximumDensity = max(self.minimumDensity, min(1.0, maximumDensity))

        if let targetVolumeFractionStart {
            self.targetVolumeFractionStart = max(self.minimumDensity, min(self.maximumDensity, targetVolumeFractionStart))
        } else {
            self.targetVolumeFractionStart = nil
        }
        if let targetVolumeFractionEnd {
            self.targetVolumeFractionEnd = max(self.minimumDensity, min(self.maximumDensity, targetVolumeFractionEnd))
        } else {
            self.targetVolumeFractionEnd = nil
        }

        self.moveLimit = max(1e-3, moveLimit)
        self.referenceElementStride = max(1, referenceElementStride)
        self.objectiveTolerance = max(1e-8, objectiveTolerance)
        self.densityChangeTolerance = max(0, densityChangeTolerance)
        self.explicitControls = explicitControls
    }
}

public struct PatchObjectiveInput2D {
    public var problem: FEMProblem2D
    public var solveResult: SolveResult2D
    public var patchElements: [Int]
    public var densities: [Float]
    public var referenceElement: Int
}

public typealias PatchObjectiveFunction2D = (PatchObjectiveInput2D) -> Float

public struct TopologyIterationResult2D {
    public var iteration: Int
    public var objective: Float
    public var averageDensity: Float
    public var targetVolumeFraction: Float
    public var totalDensityChange: Float
    public var volumeViolation: Float
    public var updatedElementCount: Int
}

public struct TopologyOptimizationResult2D {
    public var densities: [Float]
    public var densityHistory: [[Float]]
    public var history: [TopologyIterationResult2D]
    public var converged: Bool
    public var convergenceReason: String
    public var finalSolve: SolveResult2D
}
