#version 460 core

struct PointLight {
    vec3 position;
    vec3 color;
    float power;
};

in vec3 FragPos;
in vec3 Normal;
in vec2 UV;

out vec4 FragColor;

uniform vec3 ViewPos;
uniform PointLight light;
uniform sampler2D texture0;

void main()
{
    // Colore base: texture o default se UV invalidi
    vec3 baseColor;
    float eps = 0.001;
    if (abs(UV.x + 1.0) < eps && abs(UV.y + 1.0) < eps)
        baseColor = vec3(0.7, 0.7, 0.7);
    else
        baseColor = texture(texture0, UV).rgb;

    // Normalizzazione della normale
    vec3 norm = normalize(Normal);

    // Direzione della luce
    vec3 lightDir = normalize(light.position - FragPos);

    // Lambert diffuse
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * light.color * light.power;

    // Calcolo finale del colore
    vec3 color = baseColor * diffuse;

    FragColor = vec4(color, 1.0);
}
