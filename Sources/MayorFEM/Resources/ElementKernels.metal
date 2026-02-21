#include <metal_stdlib>
using namespace metal;

struct ElementGeometry {
    uint4 nodeIDs;
    float4 gradN0;
    float4 gradN1;
    float4 gradN2;
    float4 gradN3;
    float volume;
};

struct ElementState {
    float equivalentPlasticStrain;
    float damage;
};

struct MaterialConstants {
    float lambda;
    float mu;
    float yieldStress;
    float hardeningModulus;
    float damageOnset;
    float damageSlope;
    float damageCap;
    float regularization;
};

inline float3 grad3(float4 gradN) {
    return float3(gradN.x, gradN.y, gradN.z);
}

inline float3x3 outerProduct(float3 a, float3 b) {
    return float3x3(a * b.x, a * b.y, a * b.z);
}

inline float trace33(float3x3 matrix) {
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}

inline float doubleContraction(float3x3 matrix) {
    return dot(matrix[0], matrix[0]) + dot(matrix[1], matrix[1]) + dot(matrix[2], matrix[2]);
}

inline float3x3 inverseTranspose3x3(float3x3 matrix) {
    float3 a = matrix[0];
    float3 b = matrix[1];
    float3 c = matrix[2];

    float3 c0 = cross(b, c);
    float3 c1 = cross(c, a);
    float3 c2 = cross(a, b);

    float det = dot(a, c0);
    float invDet = 1.0f / max(det, 1e-12f);

    return float3x3(c0 * invDet, c1 * invDet, c2 * invDet);
}

kernel void elementResidualKernel(
    device const float3* referencePositions [[buffer(0)]],
    device const float3* displacements [[buffer(1)]],
    device const ElementGeometry* elements [[buffer(2)]],
    constant MaterialConstants& material [[buffer(3)]],
    device const ElementState* previousStates [[buffer(4)]],
    device float3* outElementNodeForces [[buffer(5)]],
    device ElementState* outTrialStates [[buffer(6)]],
    constant uint& elementCount [[buffer(7)]],
    uint id [[thread_position_in_grid]]
) {
    if (id >= elementCount) {
        return;
    }

    ElementGeometry element = elements[id];
    ElementState previousState = previousStates[id];

    uint i0 = element.nodeIDs.x;
    uint i1 = element.nodeIDs.y;
    uint i2 = element.nodeIDs.z;
    uint i3 = element.nodeIDs.w;

    float3 x0 = referencePositions[i0] + displacements[i0];
    float3 x1 = referencePositions[i1] + displacements[i1];
    float3 x2 = referencePositions[i2] + displacements[i2];
    float3 x3 = referencePositions[i3] + displacements[i3];

    float3 g0 = grad3(element.gradN0);
    float3 g1 = grad3(element.gradN1);
    float3 g2 = grad3(element.gradN2);
    float3 g3 = grad3(element.gradN3);

    float3x3 F = outerProduct(x0, g0) +
                 outerProduct(x1, g1) +
                 outerProduct(x2, g2) +
                 outerProduct(x3, g3);

    float3x3 identity = float3x3(float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1));

    float detF = determinant(F);
    if (detF < material.regularization) {
        F += (material.regularization - detF) * identity;
        detF = determinant(F);
    }
    detF = max(detF, material.regularization);

    float3x3 b = F * transpose(F);
    float3x3 tauTrial = material.mu * (b - identity) + material.lambda * log(detF) * identity;

    float meanStress = trace33(tauTrial) / 3.0f;
    float3x3 deviatoric = tauTrial - meanStress * identity;

    float qTrial = sqrt(max(1e-12f, 1.5f * doubleContraction(deviatoric)));
    float flowStress = material.yieldStress + material.hardeningModulus * previousState.equivalentPlasticStrain;

    float equivalentPlasticStrain = previousState.equivalentPlasticStrain;
    if (qTrial > flowStress) {
        float deltaGamma = (qTrial - flowStress) / (3.0f * material.mu + material.hardeningModulus + 1e-6f);
        float scale = max(0.0f, 1.0f - ((3.0f * material.mu * deltaGamma) / (qTrial + 1e-8f)));
        deviatoric *= scale;
        equivalentPlasticStrain += 0.816496580927726f * deltaGamma;
    }

    float damage = previousState.damage;
    if (equivalentPlasticStrain > material.damageOnset) {
        float threshold = max(previousState.equivalentPlasticStrain, material.damageOnset);
        float deltaEquivalentPlastic = max(0.0f, equivalentPlasticStrain - threshold);
        float deltaDamage = deltaEquivalentPlastic / max(material.damageSlope, 1e-6f);
        damage = min(material.damageCap, previousState.damage + deltaDamage);
    }

    float softening = max(0.0f, 1.0f - damage);
    float3x3 tau = softening * (deviatoric + meanStress * identity);

    float3x3 inverseTransposeF = inverseTranspose3x3(F);
    float3x3 firstPiola = tau * inverseTransposeF;

    outElementNodeForces[id * 4 + 0] = element.volume * (firstPiola * g0);
    outElementNodeForces[id * 4 + 1] = element.volume * (firstPiola * g1);
    outElementNodeForces[id * 4 + 2] = element.volume * (firstPiola * g2);
    outElementNodeForces[id * 4 + 3] = element.volume * (firstPiola * g3);

    outTrialStates[id].equivalentPlasticStrain = equivalentPlasticStrain;
    outTrialStates[id].damage = damage;
}
