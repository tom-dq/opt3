import Foundation
import simd

public enum ElementOrder2D {
    case linear
    case quadratic
}

public struct Tri6Connectivity {
    public var n0: UInt32
    public var n1: UInt32
    public var n2: UInt32
    public var n3: UInt32
    public var n4: UInt32
    public var n5: UInt32

    public init(_ n0: UInt32, _ n1: UInt32, _ n2: UInt32, _ n3: UInt32, _ n4: UInt32, _ n5: UInt32) {
        self.n0 = n0
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
        self.n5 = n5
    }

    subscript(index: Int) -> UInt32 {
        switch index {
        case 0: return n0
        case 1: return n1
        case 2: return n2
        case 3: return n3
        case 4: return n4
        case 5: return n5
        default:
            preconditionFailure("Tri6 index must be in 0...5")
        }
    }
}

public enum Element2D {
    case tri3(SIMD3<UInt32>)
    case tri6(Tri6Connectivity)
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

        let dx = length / Float(nx)
        let dy = height / Float(ny)

        var nodes: [SIMD2<Float>] = []
        nodes.reserveCapacity((nx + 1) * (ny + 1))

        for j in 0...ny {
            for i in 0...nx {
                nodes.append(SIMD2<Float>(Float(i) * dx, Float(j) * dy))
            }
        }

        func node(_ i: Int, _ j: Int) -> UInt32 {
            UInt32(j * (nx + 1) + i)
        }

        var tri3: [SIMD3<UInt32>] = []
        tri3.reserveCapacity(nx * ny * 2)

        for j in 0..<ny {
            for i in 0..<nx {
                let n00 = node(i, j)
                let n10 = node(i + 1, j)
                let n01 = node(i, j + 1)
                let n11 = node(i + 1, j + 1)

                tri3.append(SIMD3<UInt32>(n00, n10, n11))
                tri3.append(SIMD3<UInt32>(n00, n11, n01))
            }
        }

        switch order {
        case .linear:
            return Mesh2D(nodes: nodes, elements: tri3.map(Element2D.tri3))
        case .quadratic:
            return upgradeToQuadratic(nodes: nodes, tri3: tri3)
        }
    }

    public func subdivided(levels: Int) -> Mesh2D {
        if levels <= 0 {
            return self
        }

        var workingNodes = nodes
        var workingTriangles = linearizedTriangles(elements: elements)

        for _ in 0..<levels {
            var edgeMap: [EdgeKey2D: UInt32] = [:]
            var refined: [SIMD3<UInt32>] = []
            refined.reserveCapacity(workingTriangles.count * 4)

            func midpoint(_ a: UInt32, _ b: UInt32) -> UInt32 {
                let key = EdgeKey2D(a, b)
                if let id = edgeMap[key] {
                    return id
                }
                let pa = workingNodes[Int(a)]
                let pb = workingNodes[Int(b)]
                let p = 0.5 * (pa + pb)
                let id = UInt32(workingNodes.count)
                workingNodes.append(p)
                edgeMap[key] = id
                return id
            }

            for triangle in workingTriangles {
                let a = triangle.x
                let b = triangle.y
                let c = triangle.z

                let ab = midpoint(a, b)
                let bc = midpoint(b, c)
                let ca = midpoint(c, a)

                refined.append(SIMD3<UInt32>(a, ab, ca))
                refined.append(SIMD3<UInt32>(ab, b, bc))
                refined.append(SIMD3<UInt32>(ca, bc, c))
                refined.append(SIMD3<UInt32>(ab, bc, ca))
            }

            workingTriangles = refined
        }

        return Mesh2D(nodes: workingNodes, elements: workingTriangles.map(Element2D.tri3))
    }

    private static func upgradeToQuadratic(nodes: [SIMD2<Float>], tri3: [SIMD3<UInt32>]) -> Mesh2D {
        var upgradedNodes = nodes
        var edgeMidpoint: [EdgeKey2D: UInt32] = [:]
        var elements: [Element2D] = []
        elements.reserveCapacity(tri3.count)

        func midpoint(_ a: UInt32, _ b: UInt32) -> UInt32 {
            let key = EdgeKey2D(a, b)
            if let id = edgeMidpoint[key] {
                return id
            }
            let pa = upgradedNodes[Int(a)]
            let pb = upgradedNodes[Int(b)]
            let p = 0.5 * (pa + pb)
            let id = UInt32(upgradedNodes.count)
            upgradedNodes.append(p)
            edgeMidpoint[key] = id
            return id
        }

        for triangle in tri3 {
            let a = triangle.x
            let b = triangle.y
            let c = triangle.z

            let ab = midpoint(a, b)
            let bc = midpoint(b, c)
            let ca = midpoint(c, a)

            elements.append(.tri6(Tri6Connectivity(a, b, c, ab, bc, ca)))
        }

        return Mesh2D(nodes: upgradedNodes, elements: elements)
    }

    func linearizedTriangles(elements: [Element2D]) -> [SIMD3<UInt32>] {
        var triangles: [SIMD3<UInt32>] = []
        triangles.reserveCapacity(elements.count * 4)

        for element in elements {
            switch element {
            case .tri3(let tri):
                triangles.append(tri)
            case .tri6(let tri):
                triangles.append(SIMD3<UInt32>(tri.n0, tri.n3, tri.n5))
                triangles.append(SIMD3<UInt32>(tri.n3, tri.n1, tri.n4))
                triangles.append(SIMD3<UInt32>(tri.n5, tri.n4, tri.n2))
                triangles.append(SIMD3<UInt32>(tri.n3, tri.n4, tri.n5))
            }
        }

        return triangles
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
    public var prescribedDisplacements: [PrescribedDisplacement2D]
    public var controls: SolverControls

    public init(
        mesh: Mesh2D,
        material: MaterialParameters,
        thickness: Float = 1.0,
        prescribedDisplacements: [PrescribedDisplacement2D],
        controls: SolverControls = SolverControls()
    ) {
        self.mesh = mesh
        self.material = material
        self.thickness = thickness
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
    public var moveLimit: Float
    public var referenceElementStride: Int
    public var explicitControls: ExplicitSolverControls2D

    public init(
        iterations: Int = 8,
        patchRadius: Int = 1,
        minimumDensity: Float = 0.05,
        maximumDensity: Float = 1.0,
        moveLimit: Float = 0.12,
        referenceElementStride: Int = 1,
        explicitControls: ExplicitSolverControls2D = ExplicitSolverControls2D()
    ) {
        self.iterations = max(1, iterations)
        self.patchRadius = max(0, patchRadius)
        self.minimumDensity = max(0.001, min(1.0, minimumDensity))
        self.maximumDensity = max(self.minimumDensity, min(1.0, maximumDensity))
        self.moveLimit = max(1e-3, moveLimit)
        self.referenceElementStride = max(1, referenceElementStride)
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
    public var totalDensityChange: Float
}

public struct TopologyOptimizationResult2D {
    public var densities: [Float]
    public var history: [TopologyIterationResult2D]
    public var finalSolve: SolveResult2D
}
