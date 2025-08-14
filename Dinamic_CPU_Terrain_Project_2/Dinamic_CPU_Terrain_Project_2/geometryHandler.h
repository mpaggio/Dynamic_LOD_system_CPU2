#pragma once
#include "lib.h"
#include "strutture.h"

vector<float> generateSkyboxCube();
vector<float> simplePlane(int division, float width);
pair<vector<float>, vector<bool>> roadAndGrass(int division, float width, int roadWidth);
pair<vector<vec3>, vector<vec3>> generateSphericalBasesFromPositions(const vector<vec3>& basePositions);
vector<vec3> generateSphereMesh(vec3 center, float radius, int stacks = 16, int slices = 16);
pair<vector<vec3>, vector<vec3>> generateAllSpheresCPU(const vector<vec3>& positions, int tessLevel);
tuple<vector<float>, vector<vec4>, vector<float>, vector<vec4>> generatePatches(const vector<float>& plane, const vector<bool>& isRoad, int division);
void generateTessellatedPatchCPU(
    vector<Vertex>& cpuVertices,
    vector<unsigned int>& cpuIndices,
    const vec3& p0,
    const vec3& p1,
    const vec3& p2,
    const vec3& p3,
    int tessLevel,
    const TextureCPU& displacementMap,
    const vec4& edgeFlags,
    float displacementScale = 0.02f
);
tuple<vector<float>, vector<float>, vector<vec3>> generateBlocks(const vector<vec3>& positions, int subdivisions, bool isHedge);
pair<vector<float>, vector<float>> generatePatchesFromBlocks(const vector<float>& blocks, bool generateBases);
void generatePatchCPU(
    vector<Vertex>& cpuVertices,
    vector<unsigned int>& cpuIndices,
    const vec3& p0,
    const vec3& p1,
    const vec3& p2,
    const vec3& p3,
    const vec3& normal,
    const vector<vec3>& originalBlockVerts,
    const TextureCPU& displacementMap,
    float scale,
    int tessLevel
);
tuple<vector<float>, vector<vec3>> generateRoofs(const vector<vec3>& positions, const vector<float>& heights, int subdivisions);
pair<vector<float>, vector<float>> generatePatchesFromRoofs(const vector<float>& roofs, int subdivisions);
void generateRoofPatchCPU(
    vector<Vertex>& outVertices,
    vector<unsigned int>& outIndices,
    const vec3& p0,
    const vec3& p1,
    const vec3& p2,
    const vec3& p3,
    const vec3& normalFace,
    const vector<vec3>& originalVerts,
    const TextureCPU& dispTexture,
    float displacementScale = 0.02,
    int subdivision = 64
);
pair<vector<vec3>, vector<vec3>> generateLampLinesFromBases(const vector<vec3>& basePositions, const vector<vec3>& directions, vector<pair<vec3, vec3>>& verticalRods);
tuple<vector<vec3>, vector<vec3>, vector<vec3>> generateLampGeometryCPU(
    const vector<vec3>& basePositions,
    const vector<vec3>& directions,
    vector<pair<vec3, vec3>>& verticalRods,
    const vec3& referencePos,
    int circleSegments = 12,
    float radius = 0.01,
    int tessSegments = 64
);