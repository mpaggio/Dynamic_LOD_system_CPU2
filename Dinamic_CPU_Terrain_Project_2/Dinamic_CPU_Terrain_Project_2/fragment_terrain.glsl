#version 460 core

struct PointLight {
    vec3 position;
    vec3 color;
    float power;
};

uniform PointLight light;
uniform vec3 viewPos;
uniform sampler2D texture0;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;

out vec4 FragColor;

void main()
{
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light.position - FragPos);

    // Attenuazione morbida
    float distance = length(light.position - FragPos);
    float attenuation = 1.0 / (1.0 + 0.001 * distance);

    // Diffusa
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * light.color * attenuation;

    // Speculare
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = spec * light.color * attenuation;

    // Texture
    vec3 texColor = texture(texture0, TexCoord).rgb;

    // Risultato finale
    vec3 result = (diffuse + specular) * texColor;

    FragColor = vec4(result, 1.0);
}
