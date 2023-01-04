// Struct to represent a light
struct Light{
    vec3 position; // Position of the light
    vec3 color; // Color of the light
    float maxDistance; // Maximum distance the light can reach
};

// Calculates the distance from a point to the light source
float calculateLightDistance(vec3 point, Light light){
    return length(light.position - point);
}

// Calculates the color of a pixel given the intersection point, normal, and light source
vec3 calculateLighting(vec3 intersectionPoint, vec3 normal, Light light){
    // Calculate the distance from the intersection point to the light source
    float distance = calculateLightDistance(intersectionPoint, light);

    // If the distance is greater than the maximum range of the light, return black (i.e. the point is in shadow)
    if(distance > light.maxDistance) return vec3(0, 0, 0);

    // Calculate the diffuse lighting contribution
    float diffuse = max(dot(normal, normalize(light.position - intersectionPoint)), 0.0);
    vec3 diffuseContribution = light.color * diffuse;

    // Calculate the specular lighting contribution
    vec3 viewDirection = normalize(viewPosition - intersectionPoint);
    vec3 lightDirection = normalize(light.position - intersectionPoint);
    float specular = pow(max(dot(viewDirection, reflect(-lightDirection, normal)), 0.0), specularExponent);
    vec3 specularContribution = light.color * specular * specularStrength;

    // Return the sum of the diffuse and specular contributions
    return diffuseContribution + specularContribution;
}

// Returns the color of a pixel given the ray, scene objects, and light source
vec3 getColor(vec3 ray, Sphere spheres[SPHERE_COUNT], Light light){
    // ...

    // Calculate the lighting at the point of intersection
    vec3 lighting = calculateLighting(intersectionPoint, normal, light);

    // Return the color of the pixel
    return lighting + ambientColor + reflectionContribution;
}
