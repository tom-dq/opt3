import Foundation

@inline(__always)
private func dotProduct(_ lhs: [Float], _ rhs: [Float]) -> Float {
    precondition(lhs.count == rhs.count)
    var sum: Float = 0
    for index in lhs.indices {
        sum += lhs[index] * rhs[index]
    }
    return sum
}

@inline(__always)
private func l2Norm(_ vector: [Float]) -> Float {
    sqrt(max(0, dotProduct(vector, vector)))
}

@inline(__always)
private func applyDiagonalPreconditioner(_ vector: [Float], inverseDiagonal: [Float]) -> [Float] {
    precondition(vector.count == inverseDiagonal.count)
    var result = Array(repeating: Float(0), count: vector.count)
    for index in vector.indices {
        result[index] = inverseDiagonal[index] * vector[index]
    }
    return result
}

func solveSparseBiCGSTAB(
    matrix: SparseMatrix,
    rhs: [Float],
    tolerance: Float,
    maxIterations: Int
) throws -> [Float] {
    let n = rhs.count
    guard matrix.dimension == n else {
        throw FEMError.linearSolveFailure("Sparse matrix dimension mismatch.")
    }
    guard maxIterations > 0 else {
        throw FEMError.linearSolveFailure("Sparse solver requires maxIterations > 0.")
    }

    var x = Array(repeating: Float(0), count: n)
    var r = rhs
    var rHat = r

    var p = Array(repeating: Float(0), count: n)
    var v = Array(repeating: Float(0), count: n)

    var rhoPrevious: Float = 1
    var alpha: Float = 1
    var omega: Float = 1

    let rhsNorm = max(1e-12, l2Norm(rhs))
    if l2Norm(r) / rhsNorm <= tolerance {
        return x
    }

    let inverseDiagonal = matrix.diagonal.map { value in
        if abs(value) < 1e-12 {
            return Float(1)
        }
        return 1.0 / value
    }

    for iteration in 0..<maxIterations {
        _ = iteration
        var rho = dotProduct(rHat, r)
        if abs(rho) < 1e-20 {
            if l2Norm(r) / rhsNorm <= tolerance {
                return x
            }

            rHat = r
            rho = dotProduct(rHat, r)
            if abs(rho) < 1e-20 {
                throw FEMError.linearSolveFailure("BiCGSTAB breakdown (rho ~= 0).")
            }
        }

        let beta = (rho / rhoPrevious) * (alpha / omega)
        for index in 0..<n {
            p[index] = r[index] + beta * (p[index] - omega * v[index])
        }

        let y = applyDiagonalPreconditioner(p, inverseDiagonal: inverseDiagonal)
        v = matrix.multiply(y)

        let denominator = dotProduct(rHat, v)
        if abs(denominator) < 1e-20 {
            throw FEMError.linearSolveFailure("BiCGSTAB breakdown (alpha denominator ~= 0).")
        }

        alpha = rho / denominator

        var s = Array(repeating: Float(0), count: n)
        for index in 0..<n {
            s[index] = r[index] - alpha * v[index]
        }

        if l2Norm(s) / rhsNorm <= tolerance {
            for index in 0..<n {
                x[index] += alpha * y[index]
            }
            return x
        }

        let z = applyDiagonalPreconditioner(s, inverseDiagonal: inverseDiagonal)
        let t = matrix.multiply(z)

        let tDotT = dotProduct(t, t)
        if tDotT < 1e-20 {
            throw FEMError.linearSolveFailure("BiCGSTAB breakdown (t·t ~= 0).")
        }

        omega = dotProduct(t, s) / tDotT
        if abs(omega) < 1e-20 {
            throw FEMError.linearSolveFailure("BiCGSTAB breakdown (omega ~= 0).")
        }

        for index in 0..<n {
            x[index] += alpha * y[index] + omega * z[index]
            r[index] = s[index] - omega * t[index]
        }

        if l2Norm(r) / rhsNorm <= tolerance {
            return x
        }

        rhoPrevious = rho
    }

    throw FEMError.linearSolveFailure(
        "BiCGSTAB failed to converge in \(maxIterations) iterations (relative residual > \(tolerance))."
    )
}
