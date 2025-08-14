#version 460 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;

void main()
{
    vec3 transformedPos = vec3(aPos.x, aPos.z, aPos.y);
    vec3 transformedNormal = vec3(aNormal.x, aNormal.z, aNormal.y);

    FragPos = vec3(model * vec4(transformedPos, 1.0));
    Normal = normalize(mat3(transpose(inverse(model))) * transformedNormal);

    TexCoord = aTexCoord;

    gl_Position = proj * view * vec4(FragPos, 1.0);
}
