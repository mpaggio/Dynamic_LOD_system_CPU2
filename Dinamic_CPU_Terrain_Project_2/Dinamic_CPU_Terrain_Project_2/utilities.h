#pragma once
#include "lib.h"
#include "strutture.h"

float randomFloat(float min, float max);
int randomInt(int min, int max);
vec3 randomPosition(float L);
tuple<vector<vec3>, vector<vec3>, vector<vec3>, vector<vec3>> generateCityPositions(
    const vector<float>& plane,
    const vector<bool>& isRoad,
    int division,
    int numBlocks,
    int numHedges,
    int mapWidth
);
pair<vec3, vec3> getBoundingBox(const vector<float>& vertices);
pair<vec3, vec3> getBoundingBox(const vector<vec3>& vertices);
bool checkCollision(const vec3& playerPos, const vec3& minBB, const vec3& maxBB);
vec3 interpolate(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float u, float v);
vec3 calcNormal(const vec3& v0, const vec3& v1, const vec3& v2);
vec3 applyDisplacementCPU(
    const vec3& pos, 
    const vec3& normal, 
    const vec2& uv,
    const vector<unsigned char>& dispData, 
    int texWidth, 
    int texHeight, 
    float displacementScale
);
int lodFromDistance(float distance, float minDist, float maxDist, int minSubdiv, int maxSubdiv);