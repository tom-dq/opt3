import Foundation

struct SparseMatrix {
    var dimension: Int
    var rowPointers: [Int]
    var columns: [Int]
    var values: [Float]
    var diagonal: [Float]

    func multiply(_ vector: [Float]) -> [Float] {
        precondition(vector.count == dimension)
        var result = Array(repeating: Float(0), count: dimension)

        for row in 0..<dimension {
            var rowSum: Float = 0
            let start = rowPointers[row]
            let end = rowPointers[row + 1]
            if start < end {
                for index in start..<end {
                    rowSum += values[index] * vector[columns[index]]
                }
            }
            result[row] = rowSum
        }

        return result
    }

    func toDenseArray() -> [Float] {
        var dense = Array(repeating: Float(0), count: dimension * dimension)

        for row in 0..<dimension {
            let start = rowPointers[row]
            let end = rowPointers[row + 1]
            if start < end {
                for index in start..<end {
                    dense[row * dimension + columns[index]] = values[index]
                }
            }
        }

        return dense
    }
}

struct SparseMatrixBuilder {
    private let dimension: Int
    private var storage: [UInt64: Float]

    init(dimension: Int, reserve: Int = 0) {
        self.dimension = dimension
        self.storage = [:]
        self.storage.reserveCapacity(max(0, reserve))
    }

    mutating func add(row: Int, column: Int, value: Float) {
        guard abs(value) > 0 else {
            return
        }

        precondition((0..<dimension).contains(row))
        precondition((0..<dimension).contains(column))

        let key = Self.key(row: row, column: column)
        storage[key, default: 0] += value
    }

    func build() -> SparseMatrix {
        var rowEntries = Array(repeating: [(Int, Float)](), count: dimension)

        for (key, value) in storage {
            let row = Int(key >> 32)
            let column = Int(key & 0xFFFF_FFFF)
            rowEntries[row].append((column, value))
        }

        var rowPointers = Array(repeating: 0, count: dimension + 1)
        var columns: [Int] = []
        var values: [Float] = []
        columns.reserveCapacity(storage.count)
        values.reserveCapacity(storage.count)

        var diagonal = Array(repeating: Float(1), count: dimension)

        for row in 0..<dimension {
            rowPointers[row] = columns.count
            if !rowEntries[row].isEmpty {
                rowEntries[row].sort { lhs, rhs in lhs.0 < rhs.0 }
                for (column, value) in rowEntries[row] {
                    columns.append(column)
                    values.append(value)
                    if row == column, abs(value) > 1e-10 {
                        diagonal[row] = value
                    }
                }
            }
        }

        rowPointers[dimension] = columns.count

        return SparseMatrix(
            dimension: dimension,
            rowPointers: rowPointers,
            columns: columns,
            values: values,
            diagonal: diagonal
        )
    }

    private static func key(row: Int, column: Int) -> UInt64 {
        (UInt64(UInt32(row)) << 32) | UInt64(UInt32(column))
    }
}
