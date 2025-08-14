#version 460 core

layout(location = 0) in vec3 aPos;      // posizione vertice
layout(location = 1) in vec3 aNormal;   // normale vertice

out vec3 FragPos;
out vec3 Normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main() {
    FragPos = vec3(model * vec4(aPos, 1.0));         // posizione nello spazio mondo
    Normal = mat3(transpose(inverse(model))) * aNormal; // normale trasformata
    gl_Position = proj * view * vec4(FragPos, 1.0);
}
