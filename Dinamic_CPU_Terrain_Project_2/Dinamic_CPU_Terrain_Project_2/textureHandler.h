#pragma once
#include "lib.h"
#include "strutture.h"

GLuint loadSkybox();
vector<GLuint> loadAllTextures(const vector<const char*>& paths);
TextureCPU loadTextureCPU(const char* path);
GLuint loadSingleTexture(const string& path);
GLuint createFloatTexture2D(int width, int height, const vector<float>& data);
GLuint createDepthCubemapTexture();