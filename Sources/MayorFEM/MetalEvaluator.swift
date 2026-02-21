#if canImport(Metal)
import Foundation
import Metal
import simd

private struct MaterialConstants {
    var lambda: Float
    var mu: Float
    var yieldStress: Float
    var hardeningModulus: Float
    var damageOnset: Float
    var damageSlope: Float
    var damageCap: Float
    var regularization: Float
}

final class MetalElementEvaluator: ElementEvaluator {
    private let preparedMesh: PreparedMesh
    private let device: MTLDevice
    private let commandQueue: MTLCommandQueue
    private let pipelineState: MTLComputePipelineState

    private let referenceBuffer: MTLBuffer
    private let elementBuffer: MTLBuffer
    private let materialBuffer: MTLBuffer

    init?(preparedMesh: PreparedMesh, material: MaterialParameters) throws {
        guard let device = MTLCreateSystemDefaultDevice() else {
            return nil
        }
        guard let commandQueue = device.makeCommandQueue() else {
            throw FEMError.backendUnavailable("Failed to create Metal command queue.")
        }

        self.preparedMesh = preparedMesh
        self.device = device
        self.commandQueue = commandQueue

        let shaderURL = Bundle.module.url(forResource: "ElementKernels", withExtension: "metal")
        guard let shaderURL else {
            throw FEMError.backendUnavailable("Unable to load Metal shader source from bundle resources.")
        }

        let shaderSource = try String(contentsOf: shaderURL, encoding: .utf8)
        let library = try device.makeLibrary(source: shaderSource, options: nil)

        guard let function = library.makeFunction(name: "elementResidualKernel") else {
            throw FEMError.backendUnavailable("Metal shader function elementResidualKernel not found.")
        }

        self.pipelineState = try device.makeComputePipelineState(function: function)

        guard let referenceBuffer = device.makeBuffer(
            bytes: preparedMesh.nodes,
            length: MemoryLayout<SIMD3<Float>>.stride * preparedMesh.nodes.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal reference position buffer.")
        }

        guard let elementBuffer = device.makeBuffer(
            bytes: preparedMesh.elements,
            length: MemoryLayout<ElementGeometry>.stride * preparedMesh.elements.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal element buffer.")
        }

        let constants = MaterialConstants(
            lambda: material.lambda,
            mu: material.mu,
            yieldStress: material.yieldStress,
            hardeningModulus: material.hardeningModulus,
            damageOnset: material.damageOnset,
            damageSlope: material.damageSlope,
            damageCap: material.damageCap,
            regularization: 1e-4
        )

        guard let materialBuffer = device.makeBuffer(
            bytes: [constants],
            length: MemoryLayout<MaterialConstants>.stride,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal material buffer.")
        }

        self.referenceBuffer = referenceBuffer
        self.elementBuffer = elementBuffer
        self.materialBuffer = materialBuffer
    }

    func evaluate(displacements: [SIMD3<Float>], previousStates: [ElementState]) throws -> ElementEvaluation {
        let elementCount = preparedMesh.elements.count
        if displacements.count != preparedMesh.nodes.count {
            throw FEMError.invalidMesh("Displacement vector length does not match node count.")
        }
        if previousStates.count != elementCount {
            throw FEMError.invalidMesh("State vector length does not match element count.")
        }

        guard let displacementBuffer = device.makeBuffer(
            bytes: displacements,
            length: MemoryLayout<SIMD3<Float>>.stride * displacements.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal displacement buffer.")
        }

        guard let previousStateBuffer = device.makeBuffer(
            bytes: previousStates,
            length: MemoryLayout<ElementState>.stride * previousStates.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal state input buffer.")
        }

        guard let forceBuffer = device.makeBuffer(
            length: MemoryLayout<SIMD3<Float>>.stride * elementCount * 4,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal force output buffer.")
        }

        guard let trialStateBuffer = device.makeBuffer(
            length: MemoryLayout<ElementState>.stride * elementCount,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate Metal state output buffer.")
        }

        guard let commandBuffer = commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw FEMError.backendUnavailable("Failed to create Metal command buffer/encoder.")
        }

        var elementCountUInt32 = UInt32(elementCount)

        encoder.setComputePipelineState(pipelineState)
        encoder.setBuffer(referenceBuffer, offset: 0, index: 0)
        encoder.setBuffer(displacementBuffer, offset: 0, index: 1)
        encoder.setBuffer(elementBuffer, offset: 0, index: 2)
        encoder.setBuffer(materialBuffer, offset: 0, index: 3)
        encoder.setBuffer(previousStateBuffer, offset: 0, index: 4)
        encoder.setBuffer(forceBuffer, offset: 0, index: 5)
        encoder.setBuffer(trialStateBuffer, offset: 0, index: 6)
        encoder.setBytes(&elementCountUInt32, length: MemoryLayout<UInt32>.stride, index: 7)

        let executionWidth = max(1, pipelineState.threadExecutionWidth)
        let threadsPerThreadgroup = MTLSize(width: executionWidth, height: 1, depth: 1)
        let threadsPerGrid = MTLSize(width: elementCount, height: 1, depth: 1)

        encoder.dispatchThreads(threadsPerGrid, threadsPerThreadgroup: threadsPerThreadgroup)
        encoder.endEncoding()

        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()

        if commandBuffer.status != .completed {
            let message = commandBuffer.error?.localizedDescription ?? "unknown Metal command failure"
            throw FEMError.backendUnavailable("Metal kernel execution failed: \(message)")
        }

        let forcePointer = forceBuffer.contents().bindMemory(
            to: SIMD3<Float>.self,
            capacity: elementCount * 4
        )
        let forces = Array(UnsafeBufferPointer(start: forcePointer, count: elementCount * 4))

        let statePointer = trialStateBuffer.contents().bindMemory(
            to: ElementState.self,
            capacity: elementCount
        )
        let trialStates = Array(UnsafeBufferPointer(start: statePointer, count: elementCount))

        return ElementEvaluation(elementNodeForces: forces, trialStates: trialStates)
    }
}
#endif
