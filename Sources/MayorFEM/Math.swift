import Foundation
import simd

@inline(__always)
func outerProduct(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> simd_float3x3 {
    simd_float3x3(columns: (a * b.x, a * b.y, a * b.z))
}

@inline(__always)
func trace(_ matrix: simd_float3x3) -> Float {
    matrix.columns.0.x + matrix.columns.1.y + matrix.columns.2.z
}

@inline(__always)
func doubleContraction(_ matrix: simd_float3x3) -> Float {
    simd_dot(matrix.columns.0, matrix.columns.0) +
    simd_dot(matrix.columns.1, matrix.columns.1) +
    simd_dot(matrix.columns.2, matrix.columns.2)
}

func l2Norm(_ vector: [Float], indices: [Int]) -> Float {
    var sum: Float = 0
    for index in indices {
        let value = vector[index]
        sum += value * value
    }
    return sqrt(sum)
}

func dofVectorToNodalDisplacements(_ dofs: [Float], nodeCount: Int) -> [SIMD3<Float>] {
    var displacements = Array(repeating: SIMD3<Float>.zero, count: nodeCount)
    for node in 0..<nodeCount {
        let base = node * 3
        displacements[node] = SIMD3<Float>(dofs[base], dofs[base + 1], dofs[base + 2])
    }
    return displacements
}

func applyPrescribedValues(_ dofs: inout [Float], prescribed: [Int: Float]) {
    for (index, value) in prescribed {
        dofs[index] = value
    }
}
