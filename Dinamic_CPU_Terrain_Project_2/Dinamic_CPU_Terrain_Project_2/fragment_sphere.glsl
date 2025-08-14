#version 460 core

in vec3 FragPos;
in vec3 Center;

out vec4 FragColor;

struct PointLight {
    vec3 position;
    vec3 color;
    float power;
};

uniform vec3 ViewPos;
uniform PointLight light;

// Colore fisso della sfera
const vec3 baseColor = vec3(1.0, 0.8, 0.6);

void main() {
    // Normale della sfera calcolata dal centro al vertice
    vec3 normal = normalize(FragPos - Center);

    // Direzione luce
    vec3 lightDir = normalize(light.position - FragPos);

    // Diffuse
    float diff = max(dot(normal, lightDir), 0.0);

    // Specular (Phong)
    vec3 viewDir = normalize(ViewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);

    // Intensità luce modulata dal power
    vec3 ambient = 0.1 * baseColor;
    vec3 diffuse = diff * light.color * light.power * baseColor;
    vec3 specular = spec * light.color * light.power;

    FragColor = vec4(ambient + diffuse + specular, 1.0);
}
