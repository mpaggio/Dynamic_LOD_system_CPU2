#include "geometryHandler.h"
#include "utilities.h"

#define M_PI 3.14159265358979323846


// --- SKYBOX --- //
vector<float> generateSkyboxCube() {
    vector<float> skyboxVertices = vector<float>{
        -1.0f,  1.0f, -1.0f,  // fronte
        -1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,  // retro
        -1.0f, -1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

         1.0f, -1.0f, -1.0f,  // destra
         1.0f, -1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,  // sinistra
        -1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        -1.0f,  1.0f, -1.0f,  // alto
         1.0f,  1.0f, -1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f, -1.0f,  // basso
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f,  1.0f
    };

    return skyboxVertices;
}




// --- PLANE --- //
vector<float> simplePlane(int division, float width) {
	vector<float> plane;
	float triangleSide = width / division;

	for (int row = 0; row < division + 1; row++) {
		for (int col = 0; col < division + 1; col++) {
			vec3 vertex = vec3(col * triangleSide, 0.0, row * -triangleSide);
			plane.push_back(vertex.x);
			plane.push_back(vertex.z);
			plane.push_back(vertex.y);
		}
	}
	return plane;
}




// --- ROAD AND GRASS --- //
pair<vector<float>, vector<bool>> roadAndGrass(int division, float width, int roadWidth) {
    vector<float> vertices;
    vector<bool> isRoad;

    float cellSize = width / division;
    int center = division / 2; // indice centrale

    for (int row = 0; row <= division; row++) {
        for (int col = 0; col <= division; col++) {
            vec3 vertex = vec3(col * cellSize, 0.0f, row * -cellSize);

            // Definisco la croce centrale con larghezza roadWidth
            bool isRoadCell = (col >= center - roadWidth / 2 && col <= center + roadWidth / 2) ||
                (row >= center - roadWidth / 2 && row <= center + roadWidth / 2);

            vertices.push_back(vertex.x);
            vertices.push_back(vertex.z);
            vertices.push_back(vertex.y);

            isRoad.push_back(isRoadCell);
        }
    }

    return { vertices, isRoad };
}





// --- PLANE PATCHES --- //
tuple<vector<float>, vector<vec4>, vector<float>, vector<vec4>> generatePatches(const vector<float>& plane, const vector<bool>& isRoad, int division) {
    vector<float> roadPatches;
    vector<vec4> roadEdges; // x(top), y(right), z(bottom), w(left)
    vector<float> grassPatches; 
    vector<vec4> grassEdges; // x(top), y(right), z(bottom), w(left)

    int rowLength = division + 1;

    auto getPatchFlag = [&](int row, int col) -> vec4 {
        bool top = false, right = false, bottom = false, left = false;
        bool selfRoad = isRoad[row * (division + 1) + col];
        bool targetNeighbor = !selfRoad;

        // Controlla sopra (top)
        if (row > 0 && isRoad[(row - 1) * (division + 1) + col] == targetNeighbor)
            top = true;
        // Controlla destra (right)
        if (col < division - 1 && isRoad[row * (division + 1) + (col + 1)] == targetNeighbor)
            right = true;
        // Controlla sotto (bottom)
        if (row < division - 1 && isRoad[(row + 1) * (division + 1) + col] == targetNeighbor)
            bottom = true;
        // Controlla sinistra (left)
        if (col > 0 && isRoad[row * (division + 1) + (col - 1)] == targetNeighbor)
            left = true;
        
        return vec4(
            top ? 1.0f : 0.0f,
            right ? 1.0f : 0.0f,
            bottom ? 1.0f : 0.0f,
            left ? 1.0f : 0.0f
        );
    };

    for (int row = 0; row < division; ++row) {
        for (int col = 0; col < division; ++col) {
            
            int idx[4] = {
                ((row + 1) * rowLength + col) * 3,       // bottom-left
                ((row + 1) * rowLength + col + 1) * 3,   // bottom-right
                (row * rowLength + col + 1) * 3,         // top-right
                (row * rowLength + col) * 3              // top-left
            };

            bool patchIsRoad = isRoad[row * (division + 1) + col];
            vec4 flag = getPatchFlag(row, col);

            vector<vec4>& targetEdges = patchIsRoad ? roadEdges : grassEdges;
            vector<float>& targetPatches = patchIsRoad ? roadPatches : grassPatches;

            for (int i = 0; i < 4; ++i) {
                targetEdges.push_back(flag);
                targetPatches.push_back(plane[idx[i] + 0]);
                targetPatches.push_back(plane[idx[i] + 1]);
                targetPatches.push_back(plane[idx[i] + 2]);
            }
        }
    }
    return { roadPatches, roadEdges, grassPatches, grassEdges};
}




// --- CPU TESSELLATION --- //
void generateTessellatedPatchCPU(vector<Vertex>& cpuVertices, vector<unsigned int>& cpuIndices, const vec3& p0, const vec3& p1, 
    const vec3& p2, const vec3& p3, int tessLevel, const TextureCPU& displacementMap, const vec4& edgeFlags, float displacementScale) {
    
    int N = tessLevel;
    int baseIndex = cpuVertices.size();
    size_t baseIndexCount = cpuIndices.size();

    const float eps = 1.0f / float(N) * 0.75f;

    // Genera vertici con displacement
    for (int i = 0; i <= N; ++i) {
        float u = float(i) / N;
        for (int j = 0; j <= N; ++j) {
            float v = float(j) / N;

            // posizione bilineare
            vec3 pos = interpolate(p0, p1, p2, p3, u, v);

            // normale approssimata dai triangoli vicini (qui uso una normale iniziale verticale per la displacement)
            vec3 normal = vec3(0.0f, 0.0f, 1.0f);

            // UV per la displacement map
            vec2 uv(u, v);

            // Sample displacement
            float uu = clamp(uv.x, 0.0f, 1.0f);
            float vv = clamp(uv.y, 0.0f, 1.0f);
            int x = int(uu * (displacementMap.width - 1));
            int y = int(vv * (displacementMap.height - 1));
            float disp = displacementMap.data[y * displacementMap.width + x] / 255.0f;

            // Azzeramento/attenuazione ai bordi secondo i flag
            bool onLeft = (u <= eps);
            bool onRight = (u >= 1.0f - eps);
            bool onBottom = (v <= eps);
            bool onTop = (v >= 1.0f - eps);

            // se vuoi "duro" = proprio zero:
            if ((onTop && edgeFlags.x > 0.5f) ||
                (onRight && edgeFlags.y > 0.5f) ||
                (onBottom && edgeFlags.z > 0.5f) ||
                (onLeft && edgeFlags.w > 0.5f)) {
                disp = 0.0f;
            }

            // applica displacement
            pos += normal * (disp * displacementScale);

            cpuVertices.push_back({ pos, normal, uv });
        }
    }

    // Genera indici triangolari
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int idx0 = baseIndex + i * (N + 1) + j;
            int idx1 = baseIndex + (i + 1) * (N + 1) + j;
            int idx2 = baseIndex + (i + 1) * (N + 1) + (j + 1);
            int idx3 = baseIndex + i * (N + 1) + (j + 1);

            cpuIndices.push_back(idx0);
            cpuIndices.push_back(idx1);
            cpuIndices.push_back(idx2);

            cpuIndices.push_back(idx0);
            cpuIndices.push_back(idx2);
            cpuIndices.push_back(idx3);
        }
    }

    size_t indexCount = cpuIndices.size(); // Lunghezza totale dopo aver aggiunto questa patch
    for (size_t i = baseIndexCount; i < indexCount; i += 3) {
        unsigned int idx0 = cpuIndices[i];
        unsigned int idx1 = cpuIndices[i + 1];
        unsigned int idx2 = cpuIndices[i + 2];

        vec3 n = calcNormal(
            cpuVertices[idx0].pos,
            cpuVertices[idx1].pos,
            cpuVertices[idx2].pos
        );

        cpuVertices[idx0].normal += n;
        cpuVertices[idx1].normal += n;
        cpuVertices[idx2].normal += n;
    }

    // Normalizzazione delle normali della patch
    for (size_t i = baseIndex; i < cpuVertices.size(); ++i) {
        cpuVertices[i].normal = normalize(cpuVertices[i].normal);

        // Adattamento a x,z,y per il VS
        cpuVertices[i].normal = vec3(
            cpuVertices[i].normal.x,
            cpuVertices[i].normal.z,
            cpuVertices[i].normal.y
        );

        if (dot(cpuVertices[i].normal, vec3(0.0f, 1.0f, 0.0f)) < 0.0f) {
            cpuVertices[i].normal = -cpuVertices[i].normal;
        }
    }
}






// SPHERES
pair<vector<vec3>, vector<vec3>> generateSphericalBasesFromPositions(const vector<vec3>& basePositions) {
    vector<vec3> verts;
    vector<vec3> centers;

    // Generatore random e distribuzione per il raggio unico
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<float> radiusDist(0.06f, 0.1f);

    float radius = radiusDist(gen);  // scelgo un raggio casuale UNA VOLTA

    // Definizione fissa degli angoli/sfera
    vector<vec3> sphereCorners = {
        vec3(0,  1,  0),
        vec3(0,  -1,  0),
        vec3(0, 0,  1),
        vec3(0, 0,  -1),
        vec3(1,  0, 0),
        vec3(-1,  0, 0)
    };

    const int ottanteTriangles[8][3] = {
        {0, 2, 4}, // Ottante 1
        {0, 3, 4}, // Ottante 2
        {0, 3, 5}, // Ottante 3
        {0, 2, 5}, // Ottante 4
        {1, 2, 4}, // Ottante 5
        {1, 3, 4}, // Ottante 6
        {1, 3, 5}, // Ottante 7
        {1, 2, 5}  // Ottante 8
    };

    for (const vec3& topVertex : basePositions) {
        vec3 center = topVertex - vec3(0, radius, 0);  // Calcolo centro

        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 3; ++j) {
                vec3 dir = normalize(sphereCorners[ottanteTriangles[i][j]]);
                vec3 offset = dir * radius;
                verts.push_back(center + offset);
                centers.push_back(center);
            }
        }
    }

    return { verts, centers };
}


vector<vec3> generateSphereMesh(vec3 center, float radius, int stacks, int slices) {
    vector<vec3> vertices;
    vector<vec3> tempVertices;

    // Vertici grezzi
    for (int i = 0; i <= stacks; ++i) {
        float V = i / (float)stacks;
        float phi = V * pi<float>();

        for (int j = 0; j <= slices; ++j) {
            float U = j / (float)slices;
            float theta = U * two_pi<float>();

            float x = cos(theta) * sin(phi);
            float y = cos(phi);
            float z = sin(theta) * sin(phi);

            vec3 pos = vec3(x, y, z) * radius;
            tempVertices.push_back(pos);
        }
    }

    // Creazione triangoli
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            int first = i * (slices + 1) + j;
            int second = first + slices + 1;

            vertices.push_back(center + tempVertices[first]);
            vertices.push_back(center + tempVertices[second]);
            vertices.push_back(center + tempVertices[first + 1]);

            vertices.push_back(center + tempVertices[second]);
            vertices.push_back(center + tempVertices[second + 1]);
            vertices.push_back(center + tempVertices[first + 1]);
        }
    }

    return vertices;
}


pair<vector<vec3>, vector<vec3>> generateAllSpheresCPU(const vector<vec3>& positions, int tessLevel) {
    vector<vec3> allVertices;
    vector<vec3> allCenters;
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<float> radiusDist(0.06f, 0.1f);

    for (const vec3& pos : positions) {
        float radius = radiusDist(gen);
        vec3 center = pos - vec3(0, radius, 0);

        // Genera la mesh della sfera
        vector<vec3> sphereVertices = generateSphereMesh(center, radius, tessLevel, tessLevel);
        allVertices.insert(allVertices.end(), sphereVertices.begin(), sphereVertices.end());

        // Riempie i centri paralleli
        allCenters.insert(allCenters.end(), sphereVertices.size(), center);
    }

    return { allVertices, allCenters };
}



tuple<vector<float>, vector<float>, vector<vec3>> generateBlocks(const vector<vec3>& positions, int subdivisions, bool isHedge) {
    float baseSize = 1.0f;
    float minHeight = isHedge ? 0.4f : 1.0f;
    float maxHeight = isHedge ? 0.6f : 2.5f;

    vector<float> blocks;
    vector<float> heights;
    vector<vec3> baseVertices;

    blocks.reserve(positions.size() * subdivisions * 8 * 3);
    heights.reserve(positions.size());
    baseVertices.reserve(positions.size() * 8);

    srand(static_cast<unsigned int>(time(nullptr)));

    for (const vec3& baseCenter : positions) {
        float totalHeight = minHeight + static_cast<float>(rand()) / RAND_MAX * (maxHeight - minHeight);
        heights.push_back(totalHeight);

        float segmentHeight = totalHeight / subdivisions;

        float halfWidth = isHedge
            ? (0.2f + static_cast<float>(rand()) / RAND_MAX * 0.8f) / 2.0f  // max 1.0
            : baseSize / 2.0f;

        float halfDepth = isHedge
            ? (0.2f + static_cast<float>(rand()) / RAND_MAX * 0.8f) / 2.0f  // max 1.0
            : baseSize / 2.0f;

        // Calcolo vertici della base originale
        vec3 p0 = baseCenter + vec3(-halfWidth, 0.0f, -halfDepth);
        vec3 p1 = baseCenter + vec3(halfWidth, 0.0f, -halfDepth);
        vec3 p2 = baseCenter + vec3(halfWidth, 0.0f, halfDepth);
        vec3 p3 = baseCenter + vec3(-halfWidth, 0.0f, halfDepth);

        vec3 p4 = p0 + vec3(0.0f, totalHeight, 0.0f);
        vec3 p5 = p1 + vec3(0.0f, totalHeight, 0.0f);
        vec3 p6 = p2 + vec3(0.0f, totalHeight, 0.0f);
        vec3 p7 = p3 + vec3(0.0f, totalHeight, 0.0f);

        baseVertices.push_back(p0);
        baseVertices.push_back(p1);
        baseVertices.push_back(p2);
        baseVertices.push_back(p3);
        baseVertices.push_back(p4);
        baseVertices.push_back(p5);
        baseVertices.push_back(p6);
        baseVertices.push_back(p7);

        // Suddivisioni verticali
        for (int i = 0; i < subdivisions; ++i) {
            float y0 = i * segmentHeight;
            float y1 = y0 + segmentHeight;

            vec3 segmentBaseCenter = baseCenter + vec3(0.0f, y0, 0.0f);

            vec3 sp0 = segmentBaseCenter + vec3(-halfWidth, 0.0f, -halfDepth);
            vec3 sp1 = segmentBaseCenter + vec3(halfWidth, 0.0f, -halfDepth);
            vec3 sp2 = segmentBaseCenter + vec3(halfWidth, 0.0f, halfDepth);
            vec3 sp3 = segmentBaseCenter + vec3(-halfWidth, 0.0f, halfDepth);

            vec3 sp4 = sp0 + vec3(0.0f, segmentHeight, 0.0f);
            vec3 sp5 = sp1 + vec3(0.0f, segmentHeight, 0.0f);
            vec3 sp6 = sp2 + vec3(0.0f, segmentHeight, 0.0f);
            vec3 sp7 = sp3 + vec3(0.0f, segmentHeight, 0.0f);

            vec3 segmentCube[8] = { sp0, sp1, sp2, sp3, sp4, sp5, sp6, sp7 };

            for (int v = 0; v < 8; ++v) {
                blocks.push_back(segmentCube[v].x);
                blocks.push_back(segmentCube[v].y);
                blocks.push_back(segmentCube[v].z);
            }
        }
    }

    return { blocks, heights, baseVertices };
}

pair<vector<float>, vector<float>> generatePatchesFromBlocks(const vector<float>& blocks, bool generateBases) {
    vector<float> patches;
    vector<float> faceNormals;

    const int vertsPerBlock = 8;
    const int floatsPerVert = 3;
    const int floatsPerBlock = vertsPerBlock * floatsPerVert;

    for (size_t offset = 0; offset + floatsPerBlock <= blocks.size(); offset += floatsPerBlock) {

        // Estrai i vertici
        auto getVert = [&](int i) -> vec3 {
            return vec3(
                blocks[offset + i * 3 + 0],
                blocks[offset + i * 3 + 1],
                blocks[offset + i * 3 + 2]
            );
        };

        vec3 p0 = getVert(0); // -x -y -z
        vec3 p1 = getVert(1); // +x -y -z
        vec3 p2 = getVert(2); // +x -y +z
        vec3 p3 = getVert(3); // -x -y +z
        vec3 p4 = getVert(4); // -x +y -z
        vec3 p5 = getVert(5); // +x +y -z
        vec3 p6 = getVert(6); // +x +y +z
        vec3 p7 = getVert(7); // -x +y +z

        // Tutte le facce ordinate in senso orario rispetto alla normale uscente
        vec3 face0[4] = { p0, p1, p2, p3 }; // bottom (-Y)
        vec3 face1[4] = { p4, p5, p6, p7 }; // top (+Y)
        vec3 face2[4] = { p1, p5, p6, p2 }; // right (+X)
        vec3 face3[4] = { p0, p4, p7, p3 }; // left (-X)
        vec3 face4[4] = { p3, p2, p6, p7 }; // front (+Z)
        vec3 face5[4] = { p0, p1, p5, p4 }; // back (-Z)

        vec3* faces[6] = { face0, face1, face2, face3, face4, face5 };

        vec3 center = (p0 + p6) * 0.5f; // centro del cubo

        for (int f = 0; f < 6; ++f) {
            if (!generateBases) {
                if (f == 0 || f == 1) {
                    continue;
                }
            }

            vec3 v0 = faces[f][0];
            vec3 v1 = faces[f][1];
            vec3 v2 = faces[f][2];

            vec3 normal = normalize(cross(v1 - v0, v2 - v0));

            vec3 faceCenter = (v0 + v1 + v2 + faces[f][3]) * 0.25f;
            vec3 toFace = normalize(faceCenter - center);

            if (dot(normal, toFace) < 0.0f) {
                normal = -normal;
            }

            for (int v = 0; v < 4; ++v) {
                patches.push_back(faces[f][v].x);
                patches.push_back(faces[f][v].y);
                patches.push_back(faces[f][v].z);

                faceNormals.push_back(normal.x);
                faceNormals.push_back(normal.y);
                faceNormals.push_back(normal.z);
            }
        }
    }

    return { patches, faceNormals };
}


void generatePatchCPU(vector<Vertex>& cpuVertices, vector<unsigned int>& cpuIndices, const vec3& p0, 
    const vec3& p1, const vec3& p2, const vec3& p3, const vec3& normal, const vector<vec3>& originalBlockVerts, 
    const TextureCPU& displacementMap, float scale, int tessLevel){
    
    int N = tessLevel;
    int baseIndex = cpuVertices.size();

    auto remapUV = [](vec2 uv, const vec3& n) -> vec2 {
        if (abs(n.z - 1.0f) < 0.01f)        return uv;              // Front +Z
        if (abs(n.z + 1.0f) < 0.01f)        return vec2(1.0f - uv.x, uv.y); // Back -Z
        if (abs(n.x - 1.0f) < 0.01f)        return vec2(uv.y, 1.0f - uv.x); // Right +X
        if (abs(n.x + 1.0f) < 0.01f)        return vec2(1.0f - uv.y, uv.x); // Left -X
        if (abs(n.y - 1.0f) < 0.01f)        return vec2(uv.x, 1.0f - uv.y); // Top +Y
        if (abs(n.y + 1.0f) < 0.01f)        return uv;                      // Bottom -Y
        return uv;
    };

    auto sharesTwoComponents = [](const vec3& pos, const vec3& blockVert, float eps) -> bool {
        int count = 0;
        if (abs(pos.x - blockVert.x) < eps) count++;
        if (abs(pos.y - blockVert.y) < eps) count++;
        if (abs(pos.z - blockVert.z) < eps) count++;
        return count >= 2;
    };

    auto getDisplacedPos = [&](vec2 uv) -> vec3 {
        float u = uv.x;
        float v = uv.y;
        vec3 pos = (1 - u) * (1 - v) * p0 + u * (1 - v) * p1 + u * v * p2 + (1 - u) * v * p3;

        vec2 mappedUV = remapUV(uv, normal);

        float height = 0.0f;

        // disabilita displacement sui bordi del blocco
        bool onBorder = false;
        for (auto& bVert : originalBlockVerts) {
            if (sharesTwoComponents(pos, bVert, 0.001f)) {
                onBorder = true;
                break;
            }
        }
        
        if (!onBorder) {
            pos = applyDisplacementCPU(pos, normal, mappedUV,
                displacementMap.data,
                displacementMap.width,
                displacementMap.height,
                scale);
        }

        return pos;
    };

    // Genera vertici
    vector<vector<int>> idxGrid(N + 1, vector<int>(N + 1));
    for (int i = 0; i <= N; i++) {
        float u = float(i) / N;
        for (int j = 0; j <= N; j++) {
            float v = float(j) / N;
            vec2 uv(u, v);
            vec3 pos = getDisplacedPos(uv);

            cpuVertices.push_back({ pos, vec3(0.0f), remapUV(uv, normal) });
            idxGrid[i][j] = cpuVertices.size() - 1;
        }
    }

    // Genera indici triangolari
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int i0 = idxGrid[i][j];
            int i1 = idxGrid[i + 1][j];
            int i2 = idxGrid[i + 1][j + 1];
            int i3 = idxGrid[i][j + 1];

            cpuIndices.push_back(i0); cpuIndices.push_back(i1); cpuIndices.push_back(i2);
            cpuIndices.push_back(i0); cpuIndices.push_back(i2); cpuIndices.push_back(i3);
        }
    }

    // Calcolo normali stabili
    float delta = 1.0f / 64.0f;
    for (int i = 0; i <= N; i++) {
        float u = float(i) / N;
        for (int j = 0; j <= N; j++) {
            float v = float(j) / N;

            vec2 uv(u, v);
            vec3 posC = getDisplacedPos(uv);
            vec3 posT = getDisplacedPos(uv + vec2(delta, 0)); // direzione tangente
            vec3 posB = getDisplacedPos(uv + vec2(0, delta)); // direzione bitangente

            vec3 dT = posT - posC;
            vec3 dB = posB - posC;
            vec3 n = normalize(cross(dT, dB));

            if (dot(n, normal) < 0) n = -n;
            int idx = idxGrid[i][j];
            cpuVertices[idx].normal = n;
        }
    }
}




tuple<vector<float>, vector<vec3>> generateRoofs(const vector<vec3>& positions, const vector<float>& heights, int subdivisions) {
    float baseSize = 1.2f;
    float roofHeight = 0.2f;
    float shrinkFactor = 0.4f; // Riduzione della base superiore

    vector<float> roofs;
    roofs.reserve(positions.size() * subdivisions * subdivisions * 4 * 3 * 6); // riserva abbondante
    vector<vec3> baseVertices;
    baseVertices.reserve(positions.size() * 8); // 8 vertici per tetto base

    for (size_t i = 0; i < positions.size(); ++i) {
        vec3 roofBaseCenter = positions[i] + vec3(0.0f, heights[i], 0.0f);

        float halfSize = baseSize / 2.0f;
        float shrinkHalfSize = halfSize * shrinkFactor;

        // Vertici base inferiore (non suddivisi)
        vec3 lower[4] = {
            roofBaseCenter + vec3(-halfSize, 0.0f, -halfSize),
            roofBaseCenter + vec3(halfSize, 0.0f, -halfSize),
            roofBaseCenter + vec3(halfSize, 0.0f, halfSize),
            roofBaseCenter + vec3(-halfSize, 0.0f, halfSize)
        };

        // Vertici base superiore (non suddivisi)
        vec3 upper[4] = {
            roofBaseCenter + vec3(-shrinkHalfSize, roofHeight, -shrinkHalfSize),
            roofBaseCenter + vec3(shrinkHalfSize, roofHeight, -shrinkHalfSize),
            roofBaseCenter + vec3(shrinkHalfSize, roofHeight, shrinkHalfSize),
            roofBaseCenter + vec3(-shrinkHalfSize, roofHeight, shrinkHalfSize)
        };

        // Salvo i vertici originali NON suddivisi (8 vertici)
        baseVertices.push_back(lower[0]);
        baseVertices.push_back(lower[1]);
        baseVertices.push_back(lower[2]);
        baseVertices.push_back(lower[3]);
        baseVertices.push_back(upper[0]);
        baseVertices.push_back(upper[1]);
        baseVertices.push_back(upper[2]);
        baseVertices.push_back(upper[3]);

        // Tessellazione dei lati del tetto (come prima)
        for (int f = 0; f < 4; ++f) {
            for (int i_sub = 0; i_sub < subdivisions; ++i_sub) {
                for (int j_sub = 0; j_sub < subdivisions; ++j_sub) {
                    float u0 = float(i_sub) / subdivisions;
                    float v0 = float(j_sub) / subdivisions;
                    float u1 = float(i_sub + 1) / subdivisions;
                    float v1 = float(j_sub + 1) / subdivisions;

                    vec3 pA = mix(mix(lower[f], lower[(f + 1) % 4], u0), mix(upper[f], upper[(f + 1) % 4], u0), v0);
                    vec3 pB = mix(mix(lower[f], lower[(f + 1) % 4], u1), mix(upper[f], upper[(f + 1) % 4], u1), v0);
                    vec3 pC = mix(mix(lower[f], lower[(f + 1) % 4], u1), mix(upper[f], upper[(f + 1) % 4], u1), v1);
                    vec3 pD = mix(mix(lower[f], lower[(f + 1) % 4], u0), mix(upper[f], upper[(f + 1) % 4], u0), v1);

                    roofs.push_back(pA.x); roofs.push_back(pA.y); roofs.push_back(pA.z);
                    roofs.push_back(pD.x); roofs.push_back(pD.y); roofs.push_back(pD.z);
                    roofs.push_back(pC.x); roofs.push_back(pC.y); roofs.push_back(pC.z);
                    roofs.push_back(pB.x); roofs.push_back(pB.y); roofs.push_back(pB.z);
                }
            }
        }

        // Tessellazione base inferiore (come prima)
        for (int i_sub = 0; i_sub < subdivisions; ++i_sub) {
            for (int j_sub = 0; j_sub < subdivisions; ++j_sub) {
                float u0 = float(i_sub) / subdivisions;
                float v0 = float(j_sub) / subdivisions;
                float u1 = float(i_sub + 1) / subdivisions;
                float v1 = float(j_sub + 1) / subdivisions;

                vec3 p00 = mix(mix(lower[0], lower[1], u0), mix(lower[3], lower[2], u0), v0);
                vec3 p10 = mix(mix(lower[0], lower[1], u1), mix(lower[3], lower[2], u1), v0);
                vec3 p11 = mix(mix(lower[0], lower[1], u1), mix(lower[3], lower[2], u1), v1);
                vec3 p01 = mix(mix(lower[0], lower[1], u0), mix(lower[3], lower[2], u0), v1);

                roofs.insert(roofs.end(), {
                    p00.x, p00.y, p00.z,
                    p01.x, p01.y, p01.z,
                    p11.x, p11.y, p11.z,
                    p10.x, p10.y, p10.z
                });
            }
        }

        // Tessellazione base superiore (come prima)
        for (int i_sub = 0; i_sub < subdivisions; ++i_sub) {
            for (int j_sub = 0; j_sub < subdivisions; ++j_sub) {
                float u0 = float(i_sub) / subdivisions;
                float v0 = float(j_sub) / subdivisions;
                float u1 = float(i_sub + 1) / subdivisions;
                float v1 = float(j_sub + 1) / subdivisions;

                vec3 p00 = mix(mix(upper[0], upper[1], u0), mix(upper[3], upper[2], u0), v0);
                vec3 p10 = mix(mix(upper[0], upper[1], u1), mix(upper[3], upper[2], u1), v0);
                vec3 p11 = mix(mix(upper[0], upper[1], u1), mix(upper[3], upper[2], u1), v1);
                vec3 p01 = mix(mix(upper[0], upper[1], u0), mix(upper[3], upper[2], u0), v1);

                roofs.insert(roofs.end(), {
                    p00.x, p00.y, p00.z,
                    p01.x, p01.y, p01.z,
                    p11.x, p11.y, p11.z,
                    p10.x, p10.y, p10.z
                });
            }
        }
    }

    return { roofs, baseVertices };
}


pair<vector<float>, vector<float>> generatePatchesFromRoofs(const vector<float>& roofs, int subdivisions) {
    vector<float> patches;
    vector<float> faceNormals;

    const int floatsPerVert = 3;
    const int quadsPerFace = subdivisions * subdivisions;
    const int totalFaces = 6;
    const int totalQuads = totalFaces * quadsPerFace;
    const int floatsPerQuad = 4 * floatsPerVert;
    const int floatsPerRoof = totalQuads * floatsPerQuad;

    for (size_t offset = 0; offset + floatsPerRoof <= roofs.size(); offset += floatsPerRoof) {
        for (int q = 0; q < totalQuads; ++q) {
            size_t quadOffset = offset + q * floatsPerQuad;
            int faceIndex = q / quadsPerFace;

            // Estrazione dei vertici in ordine CW: p0, p3, p2, p1
            vec3 p0(roofs[quadOffset + 0], roofs[quadOffset + 1], roofs[quadOffset + 2]);
            vec3 p3(roofs[quadOffset + 9], roofs[quadOffset + 10], roofs[quadOffset + 11]);
            vec3 p2(roofs[quadOffset + 6], roofs[quadOffset + 7], roofs[quadOffset + 8]);
            vec3 p1(roofs[quadOffset + 3], roofs[quadOffset + 4], roofs[quadOffset + 5]);

            // Calcolo normale con winding CW
            vec3 normal = normalize(cross(p3 - p0, p1 - p0));

            // Flip base inferiore (deve puntare in basso)
            if (faceIndex == 4 && normal.y > 0.0f) {
                normal = -normal;
            }
            // Flip base superiore (deve puntare in alto)
            else if (faceIndex == 5 && normal.y < 0.0f) {
                normal = -normal;
            }
            // Lati inclinati: punta verso l'esterno
            else if (faceIndex >= 0 && faceIndex <= 3) {
                normal = -normal;
            }

            // Inserimento vertici in ordine CW
            patches.insert(patches.end(), {
                p0.x, p0.y, p0.z,
                p3.x, p3.y, p3.z,
                p2.x, p2.y, p2.z,
                p1.x, p1.y, p1.z
                });

            // Normali duplicate per ogni vertice
            for (int i = 0; i < 4; ++i) {
                faceNormals.insert(faceNormals.end(), { normal.x, normal.y, normal.z });
            }
        }
    }

    return { patches, faceNormals };
}


void generateRoofPatchCPU(vector<Vertex>& outVertices, vector<unsigned int>& outIndices,
    const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3,
    const vec3& normalFace, const vector<vec3>& originalVerts,
    const TextureCPU& dispTexture, float displacementScale, int subdivision) {

    float delta = 1.0f / float(subdivision);
    int baseIndex = static_cast<int>(outVertices.size());

    auto remapUV = [](vec2 uv, const vec3& n) -> vec2 {
        if (abs(n.z - 1.0f) < 0.01f)        return uv;
        if (abs(n.z + 1.0f) < 0.01f)        return vec2(1.0f - uv.x, uv.y);
        if (abs(n.x - 1.0f) < 0.01f)        return vec2(uv.y, 1.0f - uv.x);
        if (abs(n.x + 1.0f) < 0.01f)        return vec2(1.0f - uv.y, uv.x);
        if (abs(n.y - 1.0f) < 0.01f)        return vec2(uv.x, 1.0f - uv.y);
        if (abs(n.y + 1.0f) < 0.01f)        return uv;
        return uv;
    };

    auto sharesTwoComponents = [](const vec3& pos, const vec3& vert, float eps) -> bool {
        int count = 0;
        if (abs(pos.x - vert.x) < eps) count++;
        if (abs(pos.y - vert.y) < eps) count++;
        if (abs(pos.z - vert.z) < eps) count++;
        return count >= 2;
    };

    auto getDisplacedPos = [&](vec2 uv) -> vec3 {
        vec3 pos = interpolate(p0, p1, p2, p3, uv.x, uv.y);
        vec2 mappedUV = remapUV(uv, normalFace);

        // Controllo se il vertice si trova sui bordi originali
        bool onBorder = (uv.x == 0.0f || uv.x == 1.0f || uv.y == 0.0f || uv.y == 1.0f);

        if (!onBorder) {
            pos = applyDisplacementCPU(pos, normalFace, mappedUV,
                dispTexture.data, dispTexture.width, dispTexture.height, displacementScale);
        }

        return pos;
    };

    // Genera vertici con displacement condizionato
    for (int i = 0; i <= subdivision; ++i) {
        for (int j = 0; j <= subdivision; ++j) {
            float u = i / float(subdivision);
            float v = j / float(subdivision);
            vec2 uv(u, v);

            vec3 normal;
            vec3 posC = getDisplacedPos(uv);
            if (uv.x == 0.0f || uv.x == 1.0f || uv.y == 0.0f || uv.y == 1.0f) {
                normal = normalFace;
            }
            else {
                vec3 posR = getDisplacedPos(vec2(std::min(u + delta, 1.0f), v));
                vec3 posU = getDisplacedPos(vec2(u, std::min(v + delta, 1.0f)));
                normal = normalize(cross(posR - posC, posU - posC));
                if (dot(normal, normalFace) < 0) normal = -normal;
            }

            outVertices.push_back({ posC, normal, remapUV(uv, normalFace) });
        }
    }

    // Genera indici triangolari
    for (int i = 0; i < subdivision; ++i) {
        for (int j = 0; j < subdivision; ++j) {
            int idx0 = baseIndex + i * (subdivision + 1) + j;
            int idx1 = idx0 + 1;
            int idx2 = idx0 + (subdivision + 1);
            int idx3 = idx2 + 1;

            outIndices.push_back(idx0); outIndices.push_back(idx2); outIndices.push_back(idx1);
            outIndices.push_back(idx1); outIndices.push_back(idx2); outIndices.push_back(idx3);
        }
    }
}



pair<vector<vec3>, vector<vec3>> generateLampLinesFromBases(const vector<vec3>& basePositions, const vector<vec3>& directions, vector<pair<vec3, vec3>>& verticalRods) {
    vector<vec3> result;
    vector<vec3> lightPositions;

    // Altezza e larghezza costante per tutti i lampioni
    float height = 1.3f;  // valore fisso
    float width = 0.3f;
    float halfWidth = width * 0.5f;

    for (size_t i = 0; i < basePositions.size(); ++i) {
        const vec3& base = basePositions[i];
        const vec3& dir = directions[i];

        float angle = atan2(dir.x, dir.z);

        // Funzione lambda per ruotare un offset attorno a Y
        auto rotateY = [&](const vec3& offset) -> vec3 {
            return vec3(
                offset.x * cos(angle) + offset.z * sin(angle),
                offset.y,
                -offset.x * sin(angle) + offset.z * cos(angle)
            );
        };

        vec3 baseLeft = base + rotateY(vec3(-halfWidth, 0.0f, 0.0f));
        vec3 topLeft = base + rotateY(vec3(-halfWidth, height, 0.0f));
        vec3 topRight = base + rotateY(vec3(halfWidth, height, 0.0f));
        float shortLeg = height * 0.15f;
        vec3 baseRight = base + rotateY(vec3(halfWidth, height - shortLeg, 0.0f));

        // Salvo asta verticale
        verticalRods.push_back({ baseLeft, topLeft });

        // Curva 1: baseLeft - topLeft - topRight - baseRight
        result.insert(result.end(), { baseLeft, topLeft, topRight, baseRight });
        // Curva 2: topLeft - topRight - baseRight - baseRight
        result.insert(result.end(), { topLeft, topRight, baseRight, baseRight });
        // Curva 3: baseLeft - baseLeft - topLeft - topRight
        result.insert(result.end(), { baseLeft, baseLeft, topLeft, topRight });

        lightPositions.push_back(baseRight);
    }

    return { result, lightPositions };
}

static inline vec3 orthogonal(const vec3& v) {
    // vettore non parallelo per inizializzare una normale
    return (abs(v.x) > abs(v.z))
        ? normalize(vec3(-v.y, v.x, 0.0f))
        : normalize(vec3(0.0f, -v.z, v.y));
}

static inline vec3 catmullRom(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    return 0.5f * ((2.0f * p1) +
        (-p0 + p2) * t +
        (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
        (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

static inline void sampleCurveCatmullRom(const vec3 c0, const vec3 c1, const vec3 c2, const vec3 c3, int tessSegments, vector<vec3>& out) {
    out.reserve(out.size() + tessSegments + 1);
    for (int i = 0; i <= tessSegments; ++i) {
        float t = float(i) / float(tessSegments);
        out.push_back(catmullRom(c0, c1, c2, c3, t));
    }
}

static inline void addTubeAlongPolyline(
    const vector<vec3>& points,
    int circleSegments,
    float radius,
    vector<vec3>& outVerts,
    vector<vec3>& outNormals)
{
    if (points.size() < 2) return;

    const int nRings = (int)points.size();
    vector<vec3> ringNormal(nRings);
    vector<vec3> ringBinormal(nRings);

    vec3 t0 = normalize(points[1] - points[0]);
    ringNormal[0] = orthogonal(t0);
    ringBinormal[0] = normalize(cross(t0, ringNormal[0]));

    for (int i = 1; i < nRings; ++i) {
        vec3 ti = normalize(points[std::min(i + 1, nRings - 1)] - points[i - 1]);
        ringNormal[i] = orthogonal(ti);
        ringBinormal[i] = normalize(cross(ti, ringNormal[i]));
    }

    for (int i = 0; i < nRings - 1; ++i) {
        const vec3& C0 = points[i];
        const vec3& C1 = points[i + 1];

        const vec3& N0 = ringNormal[i];
        const vec3& B0 = ringBinormal[i];
        const vec3& N1 = ringNormal[i + 1];
        const vec3& B1 = ringBinormal[i + 1];

        for (int k = 0; k < circleSegments; ++k) {
            float a0 = (float)k / circleSegments * 2.0f * pi<float>();
            float a1 = (float)(k + 1) / circleSegments * 2.0f * pi<float>();

            vec3 p00 = N0 * cos(a0) + B0 * sin(a0);
            vec3 p01 = N0 * cos(a1) + B0 * sin(a1);
            vec3 p10 = N1 * cos(a0) + B1 * sin(a0);
            vec3 p11 = N1 * cos(a1) + B1 * sin(a1);

            vec3 v00 = C0 + p00 * radius;
            vec3 v01 = C0 + p01 * radius;
            vec3 v10 = C1 + p10 * radius;
            vec3 v11 = C1 + p11 * radius;

            // triangolo 1
            outVerts.push_back(v00); outVerts.push_back(v10); outVerts.push_back(v01);
            outNormals.push_back(normalize(p00));
            outNormals.push_back(normalize(p10));
            outNormals.push_back(normalize(p01));

            // triangolo 2
            outVerts.push_back(v10); outVerts.push_back(v11); outVerts.push_back(v01);
            outNormals.push_back(normalize(p10));
            outNormals.push_back(normalize(p11));
            outNormals.push_back(normalize(p01));
        }
    }
}

tuple<vector<vec3>, vector<vec3>, vector<vec3>> generateLampGeometryCPU(const vector<vec3>& basePositions, const vector<vec3>& directions,
    vector<pair<vec3, vec3>>& verticalRods, const vec3& referencePos, int circleSegments, float radius, int tessSegments) {

    vector<vec3> vertices;
    vector<vec3> normals;
    vector<vec3> lightPositions;  

    const float height = 1.3f;
    const float width = 0.3f;
    const float halfWidth = width * 0.5f;
    const float shortLeg = height * 0.15f;

    for (size_t i = 0; i < basePositions.size(); ++i) {
        const vec3& base = basePositions[i];
        const vec3& dir = directions[i];

        float minDist = 1.0f;
        float maxDist = 3.0f;
        int minTess = 4;
        int maxTess = 64;
        int minCircle = 4;
        int maxCircle = 12;
        float distance = length(referencePos - base);
        int tessSegments = lodFromDistance(distance, minDist, maxDist, minTess, maxTess);
        int circleSegments = lodFromDistance(distance, minDist, maxDist, minCircle, maxCircle);

        // orientamento attorno a Y (come negli shader)
        float angle = atan2(dir.x, dir.z);

        auto rotateY = [&](const vec3& offset) -> vec3 {
            return vec3(
                offset.x * cos(angle) + offset.z * sin(angle),
                offset.y,
                -offset.x * sin(angle) + offset.z * cos(angle)
            );
        };

        // punti chiave (stessi del tuo generatore)
        vec3 baseLeft = base + rotateY(vec3(-halfWidth, 0.0f, 0.0f));
        vec3 topLeft = base + rotateY(vec3(-halfWidth, height, 0.0f));
        vec3 topRight = base + rotateY(vec3(halfWidth, height, 0.0f));
        vec3 baseRight = base + rotateY(vec3(halfWidth, height - shortLeg, 0.0f));

        // salva lo scheletro verticale
        verticalRods.push_back({ baseLeft, topLeft });

        // 2) Curva 1: baseLeft - topLeft - topRight - baseRight
        // campiona la Catmull-Rom e estrudi tubo
        {
            vector<vec3> curve;
            sampleCurveCatmullRom(baseLeft, topLeft, topRight, baseRight, tessSegments, curve);
            addTubeAlongPolyline(curve, circleSegments, radius, vertices, normals);
        }

        // 3) Curva 2: topLeft - topRight - baseRight - baseRight
        {
            vector<vec3> curve;
            sampleCurveCatmullRom(topLeft, topRight, baseRight, baseRight, tessSegments, curve);
            addTubeAlongPolyline(curve, circleSegments, radius, vertices, normals);
        }

        // 4) Curva 3: baseLeft - baseLeft - topLeft - topRight
        {
            vector<vec3> curve;
            sampleCurveCatmullRom(baseLeft, baseLeft, topLeft, topRight, tessSegments, curve);
            addTubeAlongPolyline(curve, circleSegments, radius, vertices, normals);
        }

        // Posizione luce: ultimo punto della "curva 2", che termina in baseRight
        lightPositions.push_back(baseRight);
    }

    return { vertices, normals, lightPositions };
}
