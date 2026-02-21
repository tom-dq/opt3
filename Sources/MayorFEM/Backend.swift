import Foundation
import simd

struct ElementEvaluation {
    var elementNodeForces: [SIMD3<Float>]
    var trialStates: [ElementState]
}

protocol ElementEvaluator {
    func evaluate(displacements: [SIMD3<Float>], previousStates: [ElementState]) throws -> ElementEvaluation
}
