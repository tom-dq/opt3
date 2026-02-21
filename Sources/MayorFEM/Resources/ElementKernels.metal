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

inline float2x2 outerProduct2(float2 a, float2 b) {
    return float2x2(a * b.x, a * b.y);
}

inline float2 mul2x2(float2 c0, float2 c1, float2 v) {
    return c0 * v.x + c1 * v.y;
}

kernel void elementResidualKernel2D(
    device const float2* referencePositions [[buffer(0)]],
    device const float2* displacements [[buffer(1)]],
    device const uint* nodeIDs [[buffer(2)]],
    device const float2* gradients [[buffer(3)]],
    device const float* areas [[buffer(4)]],
    constant MaterialConstants& material [[buffer(5)]],
    constant float& thickness [[buffer(6)]],
    device const ElementState* previousStates [[buffer(7)]],
    device const float* densities [[buffer(8)]],
    constant float& densityPenalty [[buffer(9)]],
    constant float& minimumDensity [[buffer(10)]],
    device float2* outElementNodeForces [[buffer(11)]],
    device ElementState* outTrialStates [[buffer(12)]],
    device float* outVonMises [[buffer(13)]],
    device float* outStrainEnergy [[buffer(14)]],
    constant uint& elementCount [[buffer(15)]],
    uint id [[thread_position_in_grid]]
) {
    if (id >= elementCount) {
        return;
    }

    uint base = id * 3;
    uint i0 = nodeIDs[base + 0];
    uint i1 = nodeIDs[base + 1];
    uint i2 = nodeIDs[base + 2];

    float2 g0 = gradients[base + 0];
    float2 g1 = gradients[base + 1];
    float2 g2 = gradients[base + 2];

    float2 x0 = referencePositions[i0] + displacements[i0];
    float2 x1 = referencePositions[i1] + displacements[i1];
    float2 x2 = referencePositions[i2] + displacements[i2];

    float2x2 F2 = outerProduct2(x0, g0) +
                  outerProduct2(x1, g1) +
                  outerProduct2(x2, g2);

    float detF2 = determinant(F2);
    float2x2 identity2 = float2x2(float2(1.0f, 0.0f), float2(0.0f, 1.0f));
    if (detF2 < material.regularization) {
        F2 += (material.regularization - detF2) * identity2;
        detF2 = determinant(F2);
    }
    detF2 = max(detF2, material.regularization);

    float3x3 F = float3x3(
        float3(F2[0].x, F2[0].y, 0.0f),
        float3(F2[1].x, F2[1].y, 0.0f),
        float3(0.0f, 0.0f, 1.0f)
    );

    float3x3 identity = float3x3(float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1));
    float3x3 b = F * transpose(F);
    float3x3 tauTrial = material.mu * (b - identity) + material.lambda * log(detF2) * identity;

    float meanStress = trace33(tauTrial) / 3.0f;
    float3x3 deviatoric = tauTrial - meanStress * identity;

    float qTrial = sqrt(max(1e-12f, 1.5f * doubleContraction(deviatoric)));
    ElementState previousState = previousStates[id];
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

    float density = clamp(densities[id], minimumDensity, 1.0f);
    float densityScale = pow(density, max(1.0f, densityPenalty));
    float softening = max(0.0f, 1.0f - damage);
    float3x3 tau = densityScale * softening * (deviatoric + meanStress * identity);
    float3x3 deviatoricKirchhoff = densityScale * softening * deviatoric;
    float vonMises = sqrt(max(0.0f, 1.5f * doubleContraction(deviatoricKirchhoff)));

    float3x3 inverseTransposeF = inverseTranspose3x3(F);
    float3x3 firstPiola = tau * inverseTransposeF;
    float2 p0 = float2(firstPiola[0].x, firstPiola[0].y);
    float2 p1 = float2(firstPiola[1].x, firstPiola[1].y);

    float scale = thickness * areas[id];
    float2 f0 = scale * mul2x2(p0, p1, g0);
    float2 f1 = scale * mul2x2(p0, p1, g1);
    float2 f2 = scale * mul2x2(p0, p1, g2);

    outElementNodeForces[base + 0] = f0;
    outElementNodeForces[base + 1] = f1;
    outElementNodeForces[base + 2] = f2;

    outTrialStates[id].equivalentPlasticStrain =
        previousState.equivalentPlasticStrain + densityScale * (equivalentPlasticStrain - previousState.equivalentPlasticStrain);
    outTrialStates[id].damage =
        previousState.damage + densityScale * (damage - previousState.damage);

    outVonMises[id] = vonMises;

    float2 u0 = displacements[i0];
    float2 u1 = displacements[i1];
    float2 u2 = displacements[i2];
    float strainEnergy = 0.5f * (dot(f0, u0) + dot(f1, u1) + dot(f2, u2));
    outStrainEnergy[id] = fabs(strainEnergy);
}
