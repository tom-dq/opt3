import Foundation

func solveDenseLinearSystem(coefficients: [Float], rhs: [Float]) throws -> [Float] {
    let n = rhs.count
    if coefficients.count != n * n {
        throw FEMError.linearSolveFailure("Coefficient matrix shape mismatch.")
    }

    var a = coefficients.map(Double.init)
    var b = rhs.map(Double.init)

    for pivot in 0..<n {
        var maxRow = pivot
        var maxValue = abs(a[pivot * n + pivot])

        if pivot + 1 < n {
            for row in (pivot + 1)..<n {
                let value = abs(a[row * n + pivot])
                if value > maxValue {
                    maxValue = value
                    maxRow = row
                }
            }
        }

        if maxValue < 1e-14 {
            throw FEMError.linearSolveFailure("Near-singular Jacobian detected at pivot \(pivot).")
        }

        if maxRow != pivot {
            for column in pivot..<n {
                a.swapAt(pivot * n + column, maxRow * n + column)
            }
            b.swapAt(pivot, maxRow)
        }

        let pivotValue = a[pivot * n + pivot]
        if pivot + 1 < n {
            for row in (pivot + 1)..<n {
                let factor = a[row * n + pivot] / pivotValue
                if factor == 0 {
                    continue
                }
                for column in pivot..<n {
                    a[row * n + column] -= factor * a[pivot * n + column]
                }
                b[row] -= factor * b[pivot]
            }
        }
    }

    var x = Array(repeating: 0.0, count: n)
    for row in stride(from: n - 1, through: 0, by: -1) {
        var rhsValue = b[row]
        if row + 1 < n {
            for column in (row + 1)..<n {
                rhsValue -= a[row * n + column] * x[column]
            }
        }

        let diagonal = a[row * n + row]
        if abs(diagonal) < 1e-14 {
            throw FEMError.linearSolveFailure("Zero diagonal encountered during back substitution.")
        }
        x[row] = rhsValue / diagonal
    }

    return x.map(Float.init)
}
