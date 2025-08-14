#version 460 core

struct PointLight {
    vec3 position;
    vec3 color;
    float power;
};

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 ViewPos;
uniform PointLight light;

void main() {
    vec3 baseColor = vec3(0.6, 0.6, 0.6);

    // Normale normalizzata
    vec3 norm = normalize(Normal);

    // Direzione luce
    vec3 lightDir = normalize(light.position - FragPos);

    // Diffuse
    float diff = max(dot(norm, lightDir), 0.0);

    // Specular
    vec3 viewDir = normalize(ViewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);

    // Attenuazione
    float distance = length(light.position - FragPos);
    float lin = 0.002;
    float quadratic = 0.0004;
    float attenuation = 1.0 / (1.0 + lin * distance + quadratic * distance * distance);

    // Combinazione finale
    vec3 ambient = 0.1 * baseColor;
    vec3 diffuse = diff * baseColor * light.color * light.power;
    vec3 specular = spec * light.color * 0.5;

    vec3 result = (ambient + diffuse + specular) * attenuation;
    FragColor = vec4(result, 1.0);
}
