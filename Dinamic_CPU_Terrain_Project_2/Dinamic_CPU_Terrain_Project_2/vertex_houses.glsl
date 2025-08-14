#version 460 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUV;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;


out vec3 worldPos;
out vec3 normal;
out vec2 uv;

void main() {
    vec4 wp = model * vec4(aPos, 1.0);

    worldPos = wp.xyz;
    normal = normalize(mat3(transpose(inverse(model))) * aNormal);
    uv = aUV;

    gl_Position = proj * view * wp;
}
