#include "lib.h"
#include "strutture.h"
#include "utilities.h"
#include "guiHandler.h"
#include "modelLoader.h"
#include "noiseHandler.h"
#include "shaderHandler.h"
#include "cameraHandler.h"
#include "bufferHandler.h"
#include "textureHandler.h"
#include "geometryHandler.h"
#include "interactionHandler.h"

int height = 600; //Altezza della finestra
int width = 600; //Larghezza della finestra

float Theta = -90.0f; //Angolo per la rotazione orizzontale
float Phi = 0.0f; //Angolo per la rotazione verticale
float moveSpeed = 0.02;
long long startTimeMillis = 0;

bool mouseLocked = true;
bool lineMode = true;
bool mainCharacter = true;

ViewSetup SetupTelecamera;
PerspectiveSetup SetupProspettiva;

pointLight light;

extern vector<unsigned int> indices;
extern vector<BoneInfo> bone_info_walking;
extern vector<BoneInfo> bone_info_standing;
extern const aiScene* scene_walking;
extern const aiScene* scene_standing;

int main() {
    mat4 model = mat4(1.0f); //(Nessuna trasformazione) --> Qui potrei scalare, ruotare o traslare 
    mat4 view = lookAt(vec3(0.0f, 0.0f, 2.0f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f)); //(Davanti all'origine)
    mat4 proj = perspective(radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f); //(FOV: 45, ASPECT: 4.3, ZNEAR: 0.1, ZFAR: 100)
    
    //GLFW
    glfwInit(); //Inizializzazione di GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4); //Specifica a GLFW che verrà utilizzato OpenGL versione 4.x (specifica la versione maggiore)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6); //Specifica a GLFW che verrà utilizzato OpenGL versione 4.6 (specifica la versione minore)
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //Richiede un core profile di OpenGL (che esclude le funzionalità deprecate)

    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor(); // Prendi il monitor principale
    const GLFWvidmode* mode = glfwGetVideoMode(primaryMonitor); // Prendi le modalità video del monitor (risoluzione, refresh rate, ecc)
    height = mode->height;
    width = mode->width;
    GLFWwindow* window = glfwCreateWindow(width, height, "Tessellation Shader", primaryMonitor, nullptr); // Crea la finestra fullscreen con le dimensioni del monitor
    //GLFWwindow* window = glfwCreateWindow(width, height, "Tessellation Shader", nullptr, nullptr);

    if (!window) { //Gestione dell'errore
        std::cerr << "Errore nella creazione della finestra GLFW\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window); //Attiva il contesto OpenGL associato alla finestra creata, rendendo il contesto corrente per il thread in cui viene chiamata


    //CALLBACKS
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    

    //Inizializzazione di GLAD (carica i puntatori alle funzioni OpenGL)
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Errore nell'inizializzazione di GLAD\n";
        return -1;
    }


    //ILLUMINAZIONE
    light.position = { -30.0, 50.0, 30.0 };
    light.color = { 1.0,1.0,1.0 };
    light.power = 2.0f;


    //TEXTURES
    GLuint modelTexture = loadSingleTexture("Model/Knight/textures/texture_embedded_0.png");
    GLuint skyboxTexture = loadSkybox();

    TextureCPU terrainDisplacement = loadTextureCPU("Texture/Terrain/Displacement_3.png");
    GLuint terrainTexture = loadSingleTexture("Texture/Terrain/Color_3.png");
    
    TextureCPU grassDisplacement = loadTextureCPU("Texture/Grass/Displacement_3.png");
    GLuint grassTexture = loadSingleTexture("Texture/Grass/Color_3.png");
    
    TextureCPU houseDisplacement = loadTextureCPU("Texture/House/Displacement.png");
    GLuint houseTexture = loadSingleTexture("Texture/House/Color.png");
    TextureCPU houseDisplacement2 = loadTextureCPU("Texture/House/Displacement_2.png");
    GLuint houseTexture2 = loadSingleTexture("Texture/House/Color_2.png");
    
    TextureCPU roofDisplacement = loadTextureCPU("Texture/Roof/Displacement.png");
    GLuint roofTexture = loadSingleTexture("Texture/Roof/Color.png");
    TextureCPU moldDisplacement = loadTextureCPU("Texture/Grass/Displacement.png");
    GLuint moldTexture = loadSingleTexture("Texture/Grass/Color.png");


    //MODEL
    //Texture
    string path = "Model/Knight/source/castle_guard_01.fbx";
    //extractEmbeddedTextures(path, "Model/Knight/textures");
    //Walking
    path = "Model/Knight/source/Walking.fbx";
    loadModel(path, WALKING);
    ModelBufferPair walkingModelPair = INIT_MODEL_BUFFERS();
    //Standing
    path = "Model/Knight/source/Standing.fbx";
    loadModel(path, STANDING);
    ModelBufferPair standingModelPair = INIT_MODEL_BUFFERS();

    //MODEL MOVEMENT
    vec3 modelMovement = vec3(0.0f);
    vec3 previousModelMovement = vec3(0.0f);
    vec3 modelWorldPos = vec3(5.0f, 0.0f, 0.0f); //posizione assoluta del modello in world space
    float rotationAngle = 0.0f;


    //SKYBOX
    vector<float> skyboxVertices = generateSkyboxCube();
    BufferPair skyboxPair = INIT_SIMPLE_VERTEX_BUFFERS(skyboxVertices);


    //PLANE
    float map_scale = 2.0;
    int division = 12 * map_scale;
    float size = 5.0f * map_scale;
    vector<Vertex> cpuGrassVertices;
    vector<unsigned int> cpuGrassIndices;
    vector<Vertex> cpuTerrainVertices;
    vector<unsigned int> cpuTerrainIndices;
    int tessLevel = 64;
    auto planeData = roadAndGrass(division, size, 4);
    vector<float> terrainVertices = planeData.first;
    vector<bool> isRoad = planeData.second;
    auto terrainPatches = generatePatches(terrainVertices, isRoad, division);
    vector<float> roadPatches = get<0>(terrainPatches);
    vector<vec4> roadEdges = get<1>(terrainPatches);
    vector<float> grassPatches = get<2>(terrainPatches);
    vector<vec4> grassEdges = get<3>(terrainPatches);


    //BLOCKS
    int numHouses = 6;
    int numHedges = 6;
    int houseSubdivision = 2;
    vector<Vertex> cpuHouseVertices1;
    vector<unsigned int> cpuHouseIndices1;
    vector<Vertex> cpuHouseVertices2;
    vector<unsigned int> cpuHouseIndices2;
    auto positionData = generateCityPositions(terrainVertices, isRoad, division, numHouses, numHedges, size);
    vector<vec3> housePositions = get<0>(positionData);
    vector<vec3> moldPositions = get<1>(positionData);
    vector<vec3> lampPositions = get<2>(positionData);
    vector<vec3> lampDirections = get<3>(positionData);
    auto blocksData = generateBlocks(housePositions, houseSubdivision, false);
    vector<float> blocksVertices = get<0>(blocksData);
    vector<float> blocksHeights = get<1>(blocksData);
    vec3* blocksOriginalVertices = get<2>(blocksData).data();
    auto patchData = generatePatchesFromBlocks(blocksVertices, false);
    vector<float> blocksPatches = patchData.first;
    vector<float> blocksPatchesA(blocksPatches.begin(), blocksPatches.begin() + blocksPatches.size() / 2);
    vector<float> blocksPatchesB(blocksPatches.begin() + blocksPatches.size() / 2, blocksPatches.end());
    vector<float> blocksNormals = patchData.second;
    vector<float> blocksNormalsA(blocksNormals.begin(), blocksNormals.begin() + blocksNormals.size() / 2);
    vector<float> blocksNormalsB(blocksNormals.begin() + blocksNormals.size() / 2, blocksNormals.end());

    
    
    vector<pair<vec3, vec3>> bvHouses;
    int floatsPerHouse = blocksVertices.size() / numHouses;
    for (int i = 0; i < numHouses; ++i) {
        int startVertex = i * floatsPerHouse;
        int endVertex = startVertex + floatsPerHouse;

        vector<float> singleHouseVertices(blocksVertices.begin() + startVertex, blocksVertices.begin() + endVertex);
        bvHouses.push_back(getBoundingBox(singleHouseVertices));
    }
    


    //ROOFS
    vector<Vertex> cpuRoofVertices;
    vector<unsigned int> cpuRoofIndices;
    auto roofsData = generateRoofs(housePositions, blocksHeights, houseSubdivision/2);
    vector<float> roofsVertices = get<0>(roofsData);
    vec3* roofsOriginalVertices = get<1>(roofsData).data();
    auto roofsPatchData = generatePatchesFromRoofs(roofsVertices, houseSubdivision/2);
    vector<float> roofsPatches = roofsPatchData.first;
    vector<float> roofsNormals = roofsPatchData.second;
    
    vector<pair<vec3, vec3>> bvRoofs;
    int floatsPerRoof = roofsVertices.size() / numHouses;
    for (int i = 0; i < numHouses; ++i) {
        int startVertex = i * floatsPerRoof;
        int endVertex = startVertex + floatsPerRoof;

        vector<float> singleRoofVertices(roofsVertices.begin() + startVertex, roofsVertices.begin() + endVertex);
        bvRoofs.push_back(getBoundingBox(singleRoofVertices));
    }


    //MOLDS
    int moldSubdivision = 1;
    vector<Vertex> cpuMoldVertices;
    vector<unsigned int> cpuMoldIndices;
    auto moldsData = generateBlocks(moldPositions, moldSubdivision, true);
    vector<float> moldsVertices = get<0>(moldsData);
    vector<float> moldsHeights = get<1>(moldsData);
    vec3* moldsOriginalVertices = get<2>(moldsData).data();
    auto moldPatchData = generatePatchesFromBlocks(moldsVertices, true);
    vector<float> moldsPatches = moldPatchData.first;
    vector<float> moldsNormals = moldPatchData.second;

    vector<pair<vec3, vec3>> bvMolds;
    int floatsPerMold = moldsVertices.size() / numHedges;
    for (int i = 0; i < numHedges; ++i) {
        int startVertex = i * floatsPerMold;
        int endVertex = startVertex + floatsPerMold;

        vector<float> singleMoldVertices(moldsVertices.begin() + startVertex, moldsVertices.begin() + endVertex);
        bvMolds.push_back(getBoundingBox(singleMoldVertices));
    }


    //LAMP
    int numLamps = lampPositions.size();
    vector<pair<vec3, vec3>> verticalRods;
    vec3 lodReferencePos = mainCharacter ? modelWorldPos : SetupTelecamera.position;
    auto lampsData = generateLampGeometryCPU(lampPositions, lampDirections, verticalRods, lodReferencePos);
    vector<vec3> lampVertices = get<0>(lampsData);
    vector<vec3> lampNormals = get<1>(lampsData);
    vector<vec3> lightPositions = get<2>(lampsData);
    BufferPair lampPair = INIT_VEC3_BUFFERS_WITH_NORMALS(lampVertices, lampNormals);

    vector<pair<vec3, vec3>> bvLamps;
    for (auto& rod : verticalRods) {
        vector<vec3> rodVertices = { rod.first, rod.second };
        bvLamps.push_back(getBoundingBox(rodVertices));
    }


    //LAMP LIGHTS
    vector<vec3> sphereCenters;
    vector<float> sphereRadii;

    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<float> radiusDist(0.06f, 0.1f);

    for (const vec3& pos : lightPositions) {
        float radius = radiusDist(gen);
        vec3 center = pos - vec3(0, radius, 0);
        sphereCenters.push_back(center);
        sphereRadii.push_back(radius);
    }


    //SHADER PROGRAMS
    unsigned int modelProgram = createSimpleShaderProgram(
        "vertex_model.glsl",
        "fragment_model.glsl"
    );
    unsigned int skyboxProgram = createSimpleShaderProgram(
        "vertex_skybox.glsl",
        "fragment_skybox.glsl"
    );
    unsigned int terrainProgram = createSimpleShaderProgram(
        "vertex_terrain.glsl",
        "fragment_terrain.glsl"
    );
    unsigned int housesProgram = createSimpleShaderProgram(
        "vertex_houses.glsl",
        "fragment_houses.glsl"
    );
    unsigned int housesProgram2 = createSimpleShaderProgram(
        "vertex_houses.glsl",
        "fragment_houses.glsl"
    );
    unsigned int roofProgram = createSimpleShaderProgram(
        "vertex_roofs.glsl",
        "fragment_roofs.glsl"
    );
    unsigned int moldsProgram = createSimpleShaderProgram(
        "vertex_houses.glsl",
        "fragment_houses.glsl"
    );
    unsigned int lampProgram = createSimpleShaderProgram(
        "vertex_lamps.glsl",
        "fragment_lamps.glsl"
    );
    unsigned int lampLightProgram = createSimpleShaderProgram(
        "vertex_sphere.glsl",
        "fragment_sphere.glsl"
    );


    //UNIFORMS
    //Terrain program
    int modelLoc_terrain = glGetUniformLocation(terrainProgram, "model");
    int viewLoc_terrain = glGetUniformLocation(terrainProgram, "view");
    int projLoc_terrain = glGetUniformLocation(terrainProgram, "proj");
    int viewPos_terrain = glGetUniformLocation(terrainProgram, "viewPos");
    int lightPosLocTerrain = glGetUniformLocation(terrainProgram, "light.position");
    int lightColorLocTerrain = glGetUniformLocation(terrainProgram, "light.color");
    int lightPowerLocTerrain = glGetUniformLocation(terrainProgram, "light.power");
    //Houses program
    int modelLoc_houses = glGetUniformLocation(housesProgram, "model");
    int viewLoc_houses = glGetUniformLocation(housesProgram, "view");
    int projLoc_houses = glGetUniformLocation(housesProgram, "proj");
    int viewPos_houses = glGetUniformLocation(housesProgram, "ViewPos");
    int lightPosLocTerrain_houses = glGetUniformLocation(housesProgram, "light.position");
    int lightColorLocTerrain_houses = glGetUniformLocation(housesProgram, "light.color");
    int lightPowerLocTerrain_houses = glGetUniformLocation(housesProgram, "light.power");
    //Houses program 2
    int modelLoc_houses2 = glGetUniformLocation(housesProgram2, "model");
    int viewLoc_houses2 = glGetUniformLocation(housesProgram2, "view");
    int projLoc_houses2 = glGetUniformLocation(housesProgram2, "proj");
    int viewPos_houses2 = glGetUniformLocation(housesProgram2, "ViewPos");
    int lightPosLocTerrain_houses2 = glGetUniformLocation(housesProgram2, "light.position");
    int lightColorLocTerrain_houses2 = glGetUniformLocation(housesProgram2, "light.color");
    int lightPowerLocTerrain_houses2 = glGetUniformLocation(housesProgram2, "light.power");
    //Molds program
    int modelLoc_molds = glGetUniformLocation(moldsProgram, "model");
    int viewLoc_molds = glGetUniformLocation(moldsProgram, "view");
    int projLoc_molds = glGetUniformLocation(moldsProgram, "proj");
    int viewPos_molds = glGetUniformLocation(moldsProgram, "ViewPos");
    int lightPosLocTerrain_molds = glGetUniformLocation(moldsProgram, "light.position");
    int lightColorLocTerrain_molds = glGetUniformLocation(moldsProgram, "light.color");
    int lightPowerLocTerrain_molds = glGetUniformLocation(moldsProgram, "light.power");
    //Lamps program
    int model_lamps = glGetUniformLocation(lampProgram, "model");
    int view_lamps = glGetUniformLocation(lampProgram, "view");
    int proj_lamps = glGetUniformLocation(lampProgram, "proj");
    int lightPos_lamps = glGetUniformLocation(lampProgram, "light.position");
    int lightColor_lamps = glGetUniformLocation(lampProgram, "light.color");
    int lightPower_lamps = glGetUniformLocation(lampProgram, "light.power");
    int viewPos_lamps = glGetUniformLocation(lampProgram, "ViewPos");
    //Lamps light program
    int model_lampsLight = glGetUniformLocation(lampLightProgram, "model");
    int view_lampsLight = glGetUniformLocation(lampLightProgram, "view");
    int proj_lampsLight = glGetUniformLocation(lampLightProgram, "proj");
    int lightPos_lampsLight = glGetUniformLocation(lampLightProgram, "light.position");
    int lightColor_lampsLight = glGetUniformLocation(lampLightProgram, "light.color");
    int lightPower_lampsLight = glGetUniformLocation(lampLightProgram, "light.power");
    int viewPos_lampsLight = glGetUniformLocation(lampLightProgram, "ViewPos");
    //Model program
    int modelLoc = glGetUniformLocation(modelProgram, "model");
    int viewLoc = glGetUniformLocation(modelProgram, "view");
    int projLoc = glGetUniformLocation(modelProgram, "proj");
    int cameraPosLoc = glGetUniformLocation(modelProgram, "ViewPos");
    int lightPosLoc = glGetUniformLocation(modelProgram, "light.position");
    int lightColorLoc = glGetUniformLocation(modelProgram, "light.color");
    int lightPowerLoc = glGetUniformLocation(modelProgram, "light.power");
    int bonesLoc = glGetUniformLocation(modelProgram, "bones");
    //Skybox program
    int viewLocSkybox = glGetUniformLocation(skyboxProgram, "View");
    int projLocSkybox = glGetUniformLocation(skyboxProgram, "Projection");
    //Roofs program
    int modelLoc_roofs = glGetUniformLocation(roofProgram, "model");
    int viewLoc_roofs = glGetUniformLocation(roofProgram, "view");
    int projLoc_roofs = glGetUniformLocation(roofProgram, "proj");
    int lightPosLocTerrain_roofs = glGetUniformLocation(roofProgram, "light.position");
    int lightColorLocTerrain_roofs = glGetUniformLocation(roofProgram, "light.color");
    int lightPowerLocTerrain_roofs = glGetUniformLocation(roofProgram, "light.power");
    int viewPos_roofs = glGetUniformLocation(roofProgram, "ViewPos");


    //SETTINGS
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //Imposta la modalità Wireframe per vedere le suddivisioni fatte dallo shader
    glDisable(GL_CULL_FACE); //Disabilita il culling
    glEnable(GL_DEPTH_TEST); //Abilita il depth test


    //TELECAMERA
    INIT_CAMERA_PROJECTION();


    //GUI
    initializeGui(window); //Inizializza la finestra di interazione


    //TIME
    startTimeMillis = static_cast<long long>(glfwGetTime() * 1000.0);


    //MAIN LOOP
    while (!glfwWindowShouldClose(window)) {


        //SETTINGS
        long long currentTimeMillis = static_cast<long long>(glfwGetTime() * 1000.0);
        float animationTimeSec = ((float)(currentTimeMillis - startTimeMillis)) / 1000.0f;

        //CLEAR
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //LINE MODE
        if (lineMode) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //Imposta la modalità Wireframe per vedere le suddivisioni fatte dallo shader
        }
        else {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); //Imposta la modalità Fill per vedere le suddivisioni riempite
        }

        float timeValue = glfwGetTime(); //Restituisce il tempo in secondi dall'avvio
        glPatchParameteri(GL_PATCH_VERTICES, 4); //Dice a OpenGL che ogni patch ha 4 vertici





        //TERRAIN
        cpuGrassVertices.clear();
        cpuGrassIndices.clear();
        cpuTerrainVertices.clear();
        cpuTerrainIndices.clear();

        vec3 lodReferencePos = mainCharacter ? modelWorldPos : SetupTelecamera.position;
        for (size_t i = 0, p = 0; i < roadPatches.size(); i += 12, ++p) {
            vec3 p0(roadPatches[i + 0], roadPatches[i + 1], roadPatches[i + 2]);
            vec3 p1(roadPatches[i + 3], roadPatches[i + 4], roadPatches[i + 5]);
            vec3 p2(roadPatches[i + 6], roadPatches[i + 7], roadPatches[i + 8]);
            vec3 p3(roadPatches[i + 9], roadPatches[i + 10], roadPatches[i + 11]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            patchCenter = vec3(patchCenter.x, patchCenter.z, patchCenter.y);
            float dist = length(patchCenter - lodReferencePos);

            // Calcola LOD dinamico
            float minDist = 1.0f;
            float maxDist = 1.5f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 6, 44);

            vec4 edge = roadEdges[p * 4 + 0];
            generateTessellatedPatchCPU(cpuTerrainVertices, cpuTerrainIndices, p0, p1, p2, p3, tessLevel, terrainDisplacement, edge);
        }

        for (size_t i = 0, p = 0; i < grassPatches.size(); i += 12, ++p) {
            vec3 p0(grassPatches[i + 0], grassPatches[i + 1], grassPatches[i + 2]);
            vec3 p1(grassPatches[i + 3], grassPatches[i + 4], grassPatches[i + 5]);
            vec3 p2(grassPatches[i + 6], grassPatches[i + 7], grassPatches[i + 8]);
            vec3 p3(grassPatches[i + 9], grassPatches[i + 10], grassPatches[i + 11]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            patchCenter = vec3(patchCenter.x, patchCenter.z, patchCenter.y);
            float dist = length(patchCenter - lodReferencePos);

            float minDist = 1.0f;
            float maxDist = 1.5f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 6, 44);

            vec4 edge = grassEdges[p * 4 + 0];
            generateTessellatedPatchCPU(cpuGrassVertices, cpuGrassIndices, p0, p1, p2, p3, tessLevel, grassDisplacement, edge);
        }

        // Aggiorna i buffer GPU
        MeshBuffers roadBuffers = INIT_VERTEX_BUFFERS(cpuTerrainVertices, cpuTerrainIndices);
        MeshBuffers grassBuffers = INIT_VERTEX_BUFFERS(cpuGrassVertices, cpuGrassIndices);



        //BLOCKS
        cpuHouseVertices1.clear();
        cpuHouseIndices1.clear();
        cpuHouseVertices2.clear();
        cpuHouseIndices2.clear();
        for (size_t i = 0, p = 0; i < blocksPatchesA.size(); i += 12, ++p) { // 4 vertici * 3 componenti
            vec3 p0(blocksPatchesA[i + 0], blocksPatchesA[i + 1], blocksPatchesA[i + 2]);
            vec3 p1(blocksPatchesA[i + 3], blocksPatchesA[i + 4], blocksPatchesA[i + 5]);
            vec3 p2(blocksPatchesA[i + 6], blocksPatchesA[i + 7], blocksPatchesA[i + 8]);
            vec3 p3(blocksPatchesA[i + 9], blocksPatchesA[i + 10], blocksPatchesA[i + 11]);

            vec3 normal(blocksNormalsA[p * 12 + 0], blocksNormalsA[p * 12 + 1], blocksNormalsA[p * 12 + 2]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            vec2 patchCenterXZ(patchCenter.x, patchCenter.z);        // solo x e z
            vec2 refPosXZ(lodReferencePos.x, lodReferencePos.z);
            float dist = length(patchCenterXZ - refPosXZ);

            // Calcola LOD dinamico
            float minDist = 1.0f;
            float maxDist = 1.5f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 6, 44);

            int vertsPerBlock = 8;
            int segmentsPerBlock = 2; // quello che hai usato in generateBlocks
            int patchesPerSegment = 4;            // solo facce laterali
            int patchesPerBlock = patchesPerSegment * segmentsPerBlock;

            int blockIndex = p / patchesPerBlock; // questo ti dà l'indice corretto del blocco
            vector<vec3> blockVerts(get<2>(blocksData).begin() + blockIndex * vertsPerBlock,
                get<2>(blocksData).begin() + (blockIndex + 1) * vertsPerBlock);

            generatePatchCPU(cpuHouseVertices1, cpuHouseIndices1, p0, p1, p2, p3, normal, blockVerts, houseDisplacement, 0.04f, tessLevel);
        }

        for (size_t i = 0, p = 0; i < blocksPatchesB.size(); i += 12, ++p) { // 4 vertici * 3 componenti
            vec3 p0(blocksPatchesB[i + 0], blocksPatchesB[i + 1], blocksPatchesB[i + 2]);
            vec3 p1(blocksPatchesB[i + 3], blocksPatchesB[i + 4], blocksPatchesB[i + 5]);
            vec3 p2(blocksPatchesB[i + 6], blocksPatchesB[i + 7], blocksPatchesB[i + 8]);
            vec3 p3(blocksPatchesB[i + 9], blocksPatchesB[i + 10], blocksPatchesB[i + 11]);

            vec3 normal(blocksNormalsB[p * 12 + 0], blocksNormalsB[p * 12 + 1], blocksNormalsB[p * 12 + 2]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            vec2 patchCenterXZ(patchCenter.x, patchCenter.z);        // solo x e z
            vec2 refPosXZ(lodReferencePos.x, lodReferencePos.z);
            float dist = length(patchCenterXZ - refPosXZ);

            // Calcola LOD dinamico
            float minDist = 1.0f;
            float maxDist = 1.5f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 6, 44);

            int vertsPerBlock = 8;
            int segmentsPerBlock = 2; // quello che hai usato in generateBlocks
            int patchesPerSegment = 4;            // solo facce laterali
            int patchesPerBlock = patchesPerSegment * segmentsPerBlock;
            int numBlocksA = 24 / (patchesPerSegment * 2);

            int blockIndexB = (p / patchesPerBlock) + numBlocksA;
            vector<vec3> blockVertsB(get<2>(blocksData).begin() + blockIndexB * vertsPerBlock,
                get<2>(blocksData).begin() + (blockIndexB + 1) * vertsPerBlock);

            generatePatchCPU(cpuHouseVertices2, cpuHouseIndices2, p0, p1, p2, p3, normal, blockVertsB, houseDisplacement2, 0.04f, tessLevel);
        }

        MeshBuffers houseBuffers1 = INIT_VERTEX_BUFFERS(cpuHouseVertices1, cpuHouseIndices1);
        MeshBuffers houseBuffers2 = INIT_VERTEX_BUFFERS(cpuHouseVertices2, cpuHouseIndices2);



        //MOLDS
        cpuMoldVertices.clear();
        cpuMoldIndices.clear();
        for (size_t i = 0, p = 0; i < moldsPatches.size(); i += 12, ++p) { // 4 vertici * 3 componenti
            vec3 p0(moldsPatches[i + 0], moldsPatches[i + 1], moldsPatches[i + 2]);
            vec3 p1(moldsPatches[i + 3], moldsPatches[i + 4], moldsPatches[i + 5]);
            vec3 p2(moldsPatches[i + 6], moldsPatches[i + 7], moldsPatches[i + 8]);
            vec3 p3(moldsPatches[i + 9], moldsPatches[i + 10], moldsPatches[i + 11]);

            vec3 normal(moldsNormals[p * 12 + 0], moldsNormals[p * 12 + 1], moldsNormals[p * 12 + 2]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            vec2 patchCenterXZ(patchCenter.x, patchCenter.z);        // solo x e z
            vec2 refPosXZ(lodReferencePos.x, lodReferencePos.z);
            float dist = length(patchCenterXZ - refPosXZ);

            // Calcola LOD dinamico
            float minDist = 1.0f;
            float maxDist = 1.5f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 6, 44);

            int vertsPerBlock = 8;
            int patchesPerBlock = 6;            // solo facce laterali

            int blockIndex = p / patchesPerBlock; // questo ti dà l'indice corretto del blocco
            vector<vec3> moldsVerts(get<2>(moldsData).begin() + blockIndex * vertsPerBlock,
                get<2>(moldsData).begin() + (blockIndex + 1) * vertsPerBlock);

            generatePatchCPU(cpuMoldVertices, cpuMoldIndices, p0, p1, p2, p3, normal, moldsVerts, grassDisplacement, 0.02f, tessLevel);
        }
        MeshBuffers moldBuffers = INIT_VERTEX_BUFFERS(cpuMoldVertices, cpuMoldIndices);


        //ROOFS
        cpuRoofVertices.clear();
        cpuRoofIndices.clear();
        for (size_t i = 0, p = 0; i < roofsPatches.size(); i += 12, ++p) {
            vec3 p0(roofsPatches[i + 0], roofsPatches[i + 1], roofsPatches[i + 2]);
            vec3 p1(roofsPatches[i + 3], roofsPatches[i + 4], roofsPatches[i + 5]);
            vec3 p2(roofsPatches[i + 6], roofsPatches[i + 7], roofsPatches[i + 8]);
            vec3 p3(roofsPatches[i + 9], roofsPatches[i + 10], roofsPatches[i + 11]);

            vec3 normal(roofsNormals[p * 12 + 0], roofsNormals[p * 12 + 1], roofsNormals[p * 12 + 2]);

            vec3 patchCenter = (p0 + p1 + p2 + p3) * 0.25f;
            vec2 patchCenterXZ(patchCenter.x, patchCenter.z);        // solo x e z
            vec2 refPosXZ(lodReferencePos.x, lodReferencePos.z);
            float dist = length(patchCenterXZ - refPosXZ);

            // Calcola LOD dinamico
            float minDist = 0.5f;
            float maxDist = 1.0f;
            int tessLevel = lodFromDistance(dist, minDist, maxDist, 12, 64);

            generateRoofPatchCPU(cpuRoofVertices, cpuRoofIndices, p0, p1, p2, p3, normal, get<1>(roofsData), roofDisplacement, 0.02, tessLevel);

        }
        MeshBuffers roofBuffers = INIT_VERTEX_BUFFERS(cpuRoofVertices, cpuRoofIndices);


        //LAMP LIGHTS
        vector<vec3> lampLightsVertex;
        vector<vec3> lampLightsCenter;

        for (size_t i = 0; i < sphereCenters.size(); ++i) {
            float distance = length(lodReferencePos - sphereCenters[i]);

            float minDist = 1.0f;
            float maxDist = 3.0f;
            int tessLevel = lodFromDistance(distance, minDist, maxDist, 4, 16);

            // Genera la mesh con il tessLevel calcolato
            vector<vec3> sphereVerts = generateSphereMesh(sphereCenters[i], sphereRadii[i], tessLevel, tessLevel);

            // Aggiungi i vertici e i centri al buffer
            int baseIndex = static_cast<int>(lampLightsVertex.size());
            lampLightsVertex.insert(lampLightsVertex.end(), sphereVerts.begin(), sphereVerts.end());
            lampLightsCenter.insert(lampLightsCenter.end(), sphereVerts.size(), sphereCenters[i]);
        }
        BufferPair lampLightsPair = INIT_SPHERE_BUFFERS(lampLightsVertex, lampLightsCenter);


        //LAMP
        auto lampsData = generateLampGeometryCPU(lampPositions, lampDirections, verticalRods, lodReferencePos);
        vector<vec3> lampVertices = get<0>(lampsData);
        vector<vec3> lampNormals = get<1>(lampsData);
        vector<vec3> lightPositions = get<2>(lampsData);
        BufferPair lampPair = INIT_VEC3_BUFFERS_WITH_NORMALS(lampVertices, lampNormals);







        //MODEL PROGRAM
        glUseProgram(modelProgram);

        glActiveTexture(GL_TEXTURE0 + 0);
        glBindTexture(GL_TEXTURE_2D, modelTexture);
        string uniformName = "modelTexture";
        GLint location = glGetUniformLocation(modelProgram, uniformName.c_str());
        glUniform1i(location, 0);

        ModelState state;
        bool isMoving = length(modelMovement - previousModelMovement) > 0.00001f;
        state = isMoving ? WALKING : STANDING;

        //aggiornamento dell'animazione del personaggio (se presente)
        if (state == WALKING) {
            if (scene_walking && scene_walking->mNumAnimations > 0 && scene_walking->mAnimations[0]) {
                float ticksPerSecond = scene_walking->mAnimations[0]->mTicksPerSecond != 0 ? scene_walking->mAnimations[0]->mTicksPerSecond : 25.0f; //quanti tick al secondo
                float timeInTicks = animationTimeSec * ticksPerSecond; //quanti tick sono passati
                float animationTimeTicks = fmod(timeInTicks, scene_walking->mAnimations[0]->mDuration); //prendo la parte decimale dell'operazione modulo (animazione continua)
                updateBoneTransforms(animationTimeTicks, state);
            }
        }
        else {
            if (scene_standing && scene_standing->mNumAnimations > 0 && scene_standing->mAnimations[0]) {
                float ticksPerSecond = scene_standing->mAnimations[0]->mTicksPerSecond != 0 ? scene_standing->mAnimations[0]->mTicksPerSecond : 25.0f; //quanti tick al secondo
                float timeInTicks = animationTimeSec * ticksPerSecond; //quanti tick sono passati
                float animationTimeTicks = fmod(timeInTicks, scene_standing->mAnimations[0]->mDuration); //prendo la parte decimale dell'operazione modulo (animazione continua)
                updateBoneTransforms(animationTimeTicks, state);
            }
        }

        mat4 objectModel = mat4(1.0f);
        objectModel = translate(objectModel, modelWorldPos);
        objectModel = scale(objectModel, vec3(0.002f));
        objectModel = rotate(objectModel, radians(float(180)), vec3(0.0f, 1.0f, 0.0f));
        objectModel = rotate(objectModel, radians(rotationAngle), vec3(0.0f, 1.0f, 0.0f));
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, value_ptr(objectModel));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(cameraPosLoc, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLoc, 1, value_ptr(light.position));
        glUniform3fv(lightColorLoc, 1, value_ptr(light.color));
        glUniform1f(lightPowerLoc, light.power);

        mat4 boneTransforms[128];

        if (state == WALKING) {
            for (int i = 0; i < bone_info_walking.size(); i++)
                boneTransforms[i] = bone_info_walking[i].finalTransform;

            glUniformMatrix4fv(bonesLoc, bone_info_walking.size(), GL_FALSE, value_ptr(boneTransforms[0]));
        }
        else {
            for (int i = 0; i < bone_info_standing.size(); i++)
                boneTransforms[i] = bone_info_standing[i].finalTransform;

            glUniformMatrix4fv(bonesLoc, bone_info_standing.size(), GL_FALSE, value_ptr(boneTransforms[0]));
        }


        glBindVertexArray(walkingModelPair.vao);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);


        //SKYBOX
        glDepthFunc(GL_LEQUAL);       // per permettere la skybox in fondo
        glDepthMask(GL_FALSE);        // disattiva scrittura nello z-buffer

        glUseProgram(skyboxProgram);
        
        glUniform1i(glGetUniformLocation(skyboxProgram, "skybox"), 0);
        glUniformMatrix4fv(viewLocSkybox, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLocSkybox, 1, GL_FALSE, value_ptr(proj));

        glBindVertexArray(skyboxPair.vao);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);

        glDrawArrays(GL_TRIANGLES, 0, 36);
        
        glBindVertexArray(0);
        glDepthMask(GL_TRUE);         // riattiva scrittura per gli oggetti normali
        glDepthFunc(GL_LESS);         // ripristina depth test standard
 

        //TERRAIN PROGRAM (1)
        glUseProgram(terrainProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, terrainTexture);
        glUniform1i(glGetUniformLocation(terrainProgram, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_terrain, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_terrain, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_terrain, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewPos_terrain, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain, light.power);
        
        glBindVertexArray(roadBuffers.vao);
        glDrawElements(GL_TRIANGLES, cpuTerrainIndices.size(), GL_UNSIGNED_INT, 0);


        //TERRAIN PROGRAM (2)
        glUseProgram(terrainProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, grassTexture);
        glUniform1i(glGetUniformLocation(terrainProgram, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_terrain, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_terrain, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_terrain, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewPos_terrain, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain, light.power);

        glBindVertexArray(grassBuffers.vao);
        glDrawElements(GL_TRIANGLES, cpuGrassIndices.size(), GL_UNSIGNED_INT, 0);


        //HOUSES PROGRAM (1)
        glUseProgram(housesProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, houseTexture);
        glUniform1i(glGetUniformLocation(housesProgram, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_houses, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_houses, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_houses, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewPos_houses, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain_houses, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain_houses, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain_houses, light.power);
        
        glBindVertexArray(houseBuffers1.vao);
        glDrawElements(GL_TRIANGLES, cpuHouseIndices1.size(), GL_UNSIGNED_INT, 0);


        //HOUSES PROGRAM (2)
        glUseProgram(housesProgram2);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, houseTexture2);
        glUniform1i(glGetUniformLocation(housesProgram2, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_houses2, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_houses2, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_houses2, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewPos_houses2, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain_houses2, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain_houses2, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain_houses2, light.power);
        
        glBindVertexArray(houseBuffers2.vao);
        glDrawElements(GL_TRIANGLES, cpuHouseIndices2.size(), GL_UNSIGNED_INT, 0);


        //ROOFS PROGRAM
        glUseProgram(roofProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, roofTexture);
        glUniform1i(glGetUniformLocation(roofProgram, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_roofs, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_roofs, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_roofs, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewLoc_roofs, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain_roofs, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain_roofs, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain_roofs, light.power);
        
        glBindVertexArray(roofBuffers.vao);
        glDrawElements(GL_TRIANGLES, cpuRoofIndices.size(), GL_UNSIGNED_INT, 0);


        //MOLDS PROGRAM
        glUseProgram(moldsProgram);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, moldTexture);
        glUniform1i(glGetUniformLocation(moldsProgram, "texture0"), 0);

        glUniformMatrix4fv(modelLoc_molds, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(viewLoc_molds, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(projLoc_molds, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(viewPos_molds, 1, value_ptr(SetupTelecamera.position));
        glUniform3fv(lightPosLocTerrain_molds, 1, value_ptr(light.position));
        glUniform3fv(lightColorLocTerrain_molds, 1, value_ptr(light.color));
        glUniform1f(lightPowerLocTerrain_molds, light.power);
        
        glBindVertexArray(moldBuffers.vao);
        glDrawElements(GL_TRIANGLES, cpuMoldIndices.size(), GL_UNSIGNED_INT, 0);
        

        //LAMPS PROGRAM
        glUseProgram(lampProgram);

        glUniformMatrix4fv(model_lamps, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(view_lamps, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(proj_lamps, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(lightPos_lamps, 1, value_ptr(light.position));
        glUniform3fv(lightColor_lamps, 1, value_ptr(light.color));
        glUniform1f(lightPower_lamps, light.power);
        glUniform3fv(viewPos_lamps, 1, value_ptr(SetupTelecamera.position));

        glBindVertexArray(lampPair.vao);
        glDrawArrays(GL_TRIANGLES, 0, lampVertices.size());

        glPatchParameteri(GL_PATCH_VERTICES, 3);

        //LAMPS LIGHT PROGRAM
        glUseProgram(lampLightProgram);

        glUniformMatrix4fv(model_lampsLight, 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(view_lampsLight, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(proj_lampsLight, 1, GL_FALSE, value_ptr(proj));
        glUniform3fv(lightPos_lampsLight, 1, value_ptr(light.position));
        glUniform3fv(lightColor_lampsLight, 1, value_ptr(light.color));
        glUniform1f(lightPower_lampsLight, light.power);
        glUniform3fv(viewPos_lampsLight, 1, value_ptr(SetupTelecamera.position));

        glBindVertexArray(lampLightsPair.vao);
        glDrawArrays(GL_TRIANGLES, 0, lampLightsVertex.size());


        renderGui();
        glfwSwapBuffers(window); //Scambia il buffer frontale con quello posteriore
        glfwPollEvents(); //Controlla e gestisce gli eventi della finestra (input tastiera, mouse, ...)


        //MATRICI DI TRASFORMAZIONE
        view = lookAt(SetupTelecamera.position, SetupTelecamera.target, SetupTelecamera.upVector);
        proj = perspective(radians(SetupProspettiva.fovY), SetupProspettiva.aspect, SetupProspettiva.near_plane, SetupProspettiva.far_plane);
    
        auto inputResult = process_input(window);
        previousModelMovement = modelMovement;
        if (length(inputResult.first) > 0.0001f) {
            bool collision = false;
            for (auto bv : bvHouses) {
                if (checkCollision(modelWorldPos + inputResult.first, bv.first, bv.second)) {
                    collision = true;
                }
            }

            for (auto bv : bvRoofs) {
                if (checkCollision(modelWorldPos + inputResult.first, bv.first, bv.second)) {
                    collision = true;
                }
            }

            for (auto bv : bvMolds) {
                if (checkCollision(modelWorldPos + inputResult.first, bv.first, bv.second)) {
                    collision = true;
                }
            }

            for (auto bv : bvLamps) {
                if (checkCollision(modelWorldPos + inputResult.first, bv.first, bv.second)) {
                    collision = true;
                }
            }

            if (!collision) {
                modelMovement += inputResult.first;
                modelWorldPos += inputResult.first;
                rotationAngle = inputResult.second;
            }
        }
    }


    //TERMINAZIONE
    glDeleteProgram(modelProgram);
    glDeleteProgram(skyboxProgram);
    destroyGui();
    glfwDestroyWindow(window); //Elimina la finestra GLFW
    glfwTerminate(); //Termina la libreria GLFW, liberando tutte le risorse rimaste
    return 0;
}