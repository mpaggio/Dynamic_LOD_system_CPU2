#version 460 core

layout(location = 0) in vec3 aPos;      // posizione vertice
layout(location = 1) in vec3 aCenter;   // centro della sfera (per esempio utile se vuoi effetti basati sul centro)

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

out vec3 FragPos;   // posizione del frammento nello spazio mondo
out vec3 Center;    // centro della sfera per eventuali effetti

void main() {
    vec4 worldPos = model * vec4(aPos, 1.0);
    FragPos = worldPos.xyz;
    Center = aCenter; // passa il centro se serve
    gl_Position = proj * view * worldPos;
}
