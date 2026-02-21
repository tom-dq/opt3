#if canImport(Metal)
import Foundation
import Metal
import simd

private struct MaterialConstants2D {
    var lambda: Float
    var mu: Float
    var yieldStress: Float
    var hardeningModulus: Float
    var damageOnset: Float
    var damageSlope: Float
    var damageCap: Float
    var regularization: Float
}

final class MetalElementEvaluator2D: ElementEvaluator2D {
    private let preparedMesh: PreparedMesh2D
    private let device: MTLDevice
    private let commandQueue: MTLCommandQueue
    private let pipelineState: MTLComputePipelineState

    private let referenceBuffer: MTLBuffer
    private let nodeIDBuffer: MTLBuffer
    private let gradientBuffer: MTLBuffer
    private let areaBuffer: MTLBuffer
    private let materialBuffer: MTLBuffer
    private let thickness: Float

    init?(preparedMesh: PreparedMesh2D, material: MaterialParameters, thickness: Float) throws {
        guard let device = MTLCreateSystemDefaultDevice() else {
            return nil
        }
        guard let commandQueue = device.makeCommandQueue() else {
            throw FEMError.backendUnavailable("Failed to create Metal command queue for 2D evaluator.")
        }

        self.preparedMesh = preparedMesh
        self.device = device
        self.commandQueue = commandQueue
        self.thickness = thickness

        let shaderURL = Bundle.module.url(forResource: "ElementKernels", withExtension: "metal")
        guard let shaderURL else {
            throw FEMError.backendUnavailable("Unable to load Metal shader source from bundle resources.")
        }

        let shaderSource = try String(contentsOf: shaderURL, encoding: .utf8)
        let library = try device.makeLibrary(source: shaderSource, options: nil)

        guard let function = library.makeFunction(name: "elementResidualKernel2D") else {
            throw FEMError.backendUnavailable("Metal shader function elementResidualKernel2D not found.")
        }

        self.pipelineState = try device.makeComputePipelineState(function: function)

        var nodeIDs: [UInt32] = []
        var gradients: [SIMD2<Float>] = []
        var areas: [Float] = []
        nodeIDs.reserveCapacity(preparedMesh.elements.count * 3)
        gradients.reserveCapacity(preparedMesh.elements.count * 3)
        areas.reserveCapacity(preparedMesh.elements.count)

        for element in preparedMesh.elements {
            nodeIDs.append(element.nodeIDs[0])
            nodeIDs.append(element.nodeIDs[1])
            nodeIDs.append(element.nodeIDs[2])
            gradients.append(element.gradN0)
            gradients.append(element.gradN1)
            gradients.append(element.gradN2)
            areas.append(element.area)
        }

        guard let referenceBuffer = device.makeBuffer(
            bytes: preparedMesh.nodes,
            length: MemoryLayout<SIMD2<Float>>.stride * preparedMesh.nodes.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D reference position buffer.")
        }

        guard let nodeIDBuffer = device.makeBuffer(
            bytes: nodeIDs,
            length: MemoryLayout<UInt32>.stride * nodeIDs.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D node-ID buffer.")
        }

        guard let gradientBuffer = device.makeBuffer(
            bytes: gradients,
            length: MemoryLayout<SIMD2<Float>>.stride * gradients.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D gradient buffer.")
        }

        guard let areaBuffer = device.makeBuffer(
            bytes: areas,
            length: MemoryLayout<Float>.stride * areas.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D area buffer.")
        }

        let constants = MaterialConstants2D(
            lambda: material.lambda,
            mu: material.mu,
            yieldStress: material.yieldStress,
            hardeningModulus: material.hardeningModulus,
            damageOnset: material.damageOnset,
            damageSlope: material.damageSlope,
            damageCap: material.damageCap,
            regularization: 1e-5
        )

        guard let materialBuffer = device.makeBuffer(
            bytes: [constants],
            length: MemoryLayout<MaterialConstants2D>.stride,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D material buffer.")
        }

        self.referenceBuffer = referenceBuffer
        self.nodeIDBuffer = nodeIDBuffer
        self.gradientBuffer = gradientBuffer
        self.areaBuffer = areaBuffer
        self.materialBuffer = materialBuffer
    }

    func evaluate(
        displacements: [SIMD2<Float>],
        previousStates: [ElementState],
        densities: [Float],
        densityPenalty: Float,
        minimumDensity: Float
    ) throws -> ElementEvaluation2D {
        let elementCount = preparedMesh.elements.count

        if displacements.count != preparedMesh.nodes.count {
            throw FEMError.invalidMesh("2D displacement vector length does not match node count.")
        }
        if previousStates.count != elementCount {
            throw FEMError.invalidMesh("2D state vector length does not match element count.")
        }
        if densities.count != elementCount {
            throw FEMError.invalidMesh("2D density vector length does not match element count.")
        }

        guard let displacementBuffer = device.makeBuffer(
            bytes: displacements,
            length: MemoryLayout<SIMD2<Float>>.stride * displacements.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D displacement buffer.")
        }

        guard let previousStateBuffer = device.makeBuffer(
            bytes: previousStates,
            length: MemoryLayout<ElementState>.stride * previousStates.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D previous-state buffer.")
        }

        guard let densityBuffer = device.makeBuffer(
            bytes: densities,
            length: MemoryLayout<Float>.stride * densities.count,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D density buffer.")
        }

        guard let forceBuffer = device.makeBuffer(
            length: MemoryLayout<SIMD2<Float>>.stride * elementCount * 3,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D force output buffer.")
        }

        guard let trialStateBuffer = device.makeBuffer(
            length: MemoryLayout<ElementState>.stride * elementCount,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D state output buffer.")
        }

        guard let vonMisesBuffer = device.makeBuffer(
            length: MemoryLayout<Float>.stride * elementCount,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D von-Mises output buffer.")
        }

        guard let strainEnergyBuffer = device.makeBuffer(
            length: MemoryLayout<Float>.stride * elementCount,
            options: .storageModeShared
        ) else {
            throw FEMError.backendUnavailable("Failed to allocate 2D strain-energy output buffer.")
        }

        guard let commandBuffer = commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw FEMError.backendUnavailable("Failed to create 2D Metal command buffer/encoder.")
        }

        var elementCountUInt32 = UInt32(elementCount)
        var localThickness = thickness
        var localDensityPenalty = densityPenalty
        var localMinimumDensity = minimumDensity

        encoder.setComputePipelineState(pipelineState)
        encoder.setBuffer(referenceBuffer, offset: 0, index: 0)
        encoder.setBuffer(displacementBuffer, offset: 0, index: 1)
        encoder.setBuffer(nodeIDBuffer, offset: 0, index: 2)
        encoder.setBuffer(gradientBuffer, offset: 0, index: 3)
        encoder.setBuffer(areaBuffer, offset: 0, index: 4)
        encoder.setBuffer(materialBuffer, offset: 0, index: 5)
        encoder.setBytes(&localThickness, length: MemoryLayout<Float>.stride, index: 6)
        encoder.setBuffer(previousStateBuffer, offset: 0, index: 7)
        encoder.setBuffer(densityBuffer, offset: 0, index: 8)
        encoder.setBytes(&localDensityPenalty, length: MemoryLayout<Float>.stride, index: 9)
        encoder.setBytes(&localMinimumDensity, length: MemoryLayout<Float>.stride, index: 10)
        encoder.setBuffer(forceBuffer, offset: 0, index: 11)
        encoder.setBuffer(trialStateBuffer, offset: 0, index: 12)
        encoder.setBuffer(vonMisesBuffer, offset: 0, index: 13)
        encoder.setBuffer(strainEnergyBuffer, offset: 0, index: 14)
        encoder.setBytes(&elementCountUInt32, length: MemoryLayout<UInt32>.stride, index: 15)

        let executionWidth = max(1, pipelineState.threadExecutionWidth)
        let threadsPerThreadgroup = MTLSize(width: executionWidth, height: 1, depth: 1)
        let threadsPerGrid = MTLSize(width: elementCount, height: 1, depth: 1)

        encoder.dispatchThreads(threadsPerGrid, threadsPerThreadgroup: threadsPerThreadgroup)
        encoder.endEncoding()

        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()

        if commandBuffer.status != .completed {
            let message = commandBuffer.error?.localizedDescription ?? "unknown Metal command failure"
            throw FEMError.backendUnavailable("2D Metal kernel execution failed: \(message)")
        }

        let forcePointer = forceBuffer.contents().bindMemory(
            to: SIMD2<Float>.self,
            capacity: elementCount * 3
        )
        let forces = Array(UnsafeBufferPointer(start: forcePointer, count: elementCount * 3))

        let statePointer = trialStateBuffer.contents().bindMemory(
            to: ElementState.self,
            capacity: elementCount
        )
        let trialStates = Array(UnsafeBufferPointer(start: statePointer, count: elementCount))

        let vonMisesPointer = vonMisesBuffer.contents().bindMemory(
            to: Float.self,
            capacity: elementCount
        )
        let vonMises = Array(UnsafeBufferPointer(start: vonMisesPointer, count: elementCount))

        let energyPointer = strainEnergyBuffer.contents().bindMemory(
            to: Float.self,
            capacity: elementCount
        )
        let energies = Array(UnsafeBufferPointer(start: energyPointer, count: elementCount))

        return ElementEvaluation2D(
            elementNodeForces: forces,
            trialStates: trialStates,
            elementVonMises: vonMises,
            elementEnergies: energies
        )
    }
}
#endif
