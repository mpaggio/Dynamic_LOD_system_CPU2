#version 460 core

in vec3 worldPos;
in vec3 normal;
in vec2 uv;

out vec4 FragColor;

struct PointLight {
    vec3 position;
    vec3 color;
    float power;
};

uniform vec3 ViewPos;
uniform PointLight light;
uniform sampler2D texture0; // color

void main() {
    vec3 N = normalize(normal);
    vec3 L = normalize(light.position - worldPos);
    vec3 V = normalize(ViewPos - worldPos);
    vec3 R = reflect(-L, N);

    float ambientStrength = 0.1;
    float specularStrength = 0.5;

    vec3 ambient = ambientStrength * light.color;
    float diff = max(dot(N, L), 0.0);
    vec3 diffuse = diff * light.color;
    float spec = pow(max(dot(V, R), 0.0), 32.0);
    vec3 specular = specularStrength * spec * light.color;

    vec3 baseColor = texture(texture0, uv).rgb;
    vec3 lighting = (ambient + diffuse + specular) * baseColor * light.power;

    FragColor = vec4(lighting, 1.0);
}
