#define ANTIALIASING

// Constants for the scene configuration
const int SPHERE_COUNT = 10; // Number of spheres to render
const int DEPTH = 3; // Depth of recursion for reflections
const float START_FOV = 120.0; // Starting field of view for the camera
const float VARIANCE_FOV = 20.0; // Variance in field of view for the camera
const float CAMERA_HEIGHT = 2.0; // Height of the camera above the floor

const float WALL_REFLECTIVENESS = 0.1; // Reflectiveness of the walls
const float CEILING_REFLECTIVENESS = 0.0; // Reflectiveness of the walls
const float FLOOR_REFLECTIVENESS = 0.01; // Reflectiveness of the walls
const float SPHERE_REFLECTIVENESS = 0.3; // Reflectiveness of the spheres

vec3 backgroundColor = vec3(0.2,0.2,0.2); // Background color

// Calculates the distance from the projection plane (the plane that the camera's rays
// will be traced through) to the center of the camera's view
float calculateDistanceFromProjectionPlane(float windowRatio, float fov){
    // Calculate the tangent of half the field of view
    float tangent = 1.0/windowRatio;
    float alpha = tan(radians(fov/2.0));
    // Calculate the distance from the projection plane
    float distanceFromPlane = tangent/alpha;
    // Return the absolute value of the distance (since it could be negative)
    return abs(distanceFromPlane);
}

// Struct to represent the camera
struct Camera{
    vec3 position; // Position of the camera
    vec2 rotation; // Rotation of the camera (in x and y)
    float windowRatio; // Ratio of the width to the height of the output window
    float distanceFromProjectionPlane; // Distance from the projection plane to the center of the camera's view
};

// Struct to represent a sphere
struct Sphere{
    vec3 position; // Position of the center of the sphere
    float radius; // Radius of the sphere
    vec3 color; // Color of the sphere
    float ambientStrength; // Strength of the ambient lighting on the sphere
    float specularStrength; // Strength of the specular lighting on the sphere
    float specularExponent; // Exponent for the specular lighting on the sphere
    float reflectiveness; // Reflectiveness of the sphere
};

// Struct to represent a light
struct Light{
    vec3 position; // Position of the light
    vec3 color; // Color of the light
    float maxDistance; // Maximum distance the light can reach
};

// ==================================================================

// Calculates the distance from a point to the light source
float calculateLightDistance(vec3 point, Light light){
    return length(light.position - point);
}

// Calculates the color of a pixel given the intersection point, normal, and light source
vec3 calculateLighting(vec3 intersectionPoint, vec3 normal, Light light, float specularExponent, float specularStrength, Camera camera){
    // Calculate the distance from the intersection point to the light source
    float distance = calculateLightDistance(intersectionPoint, light);

    // If the distance is greater than the maximum range of the light, return black (i.e. the point is in shadow)
    if(distance > light.maxDistance) return vec3(0, 0, 0);

    // Calculate the diffuse lighting contribution
    float diffuse = max(dot(normal, normalize(light.position - intersectionPoint)), 0.0);
    vec3 diffuseContribution = light.color * diffuse;

    // Calculate the specular lighting contribution
    vec3 viewDirection = normalize(camera.position - intersectionPoint);
    vec3 lightDirection = normalize(light.position - intersectionPoint);
    float specular = pow(max(dot(viewDirection, reflect(-lightDirection, normal)), 0.0), specularExponent);
    vec3 specularContribution = light.color * specular * specularStrength;

    // Return the sum of the diffuse and specular contributions
    return diffuseContribution + specularContribution;
}

// Returns the color of a pixel given the ray, scene objects, and light source
vec3 getColor(vec3 ray, Sphere spheres[SPHERE_COUNT], Light light, Camera camera, vec3 ambientColor, vec3 reflectionContribution, vec3 intersectionPoint, vec3 normal, float specularExponent, float specularStrength){
    // Calculate the lighting at the point of intersection
    vec3 lighting = calculateLighting(intersectionPoint, normal, light, specularExponent, specularStrength, camera);

    // Return the color of the pixel
    return lighting + ambientColor + reflectionContribution;
}


// ==================================================================

// Returns a 3D vector representing the direction of the ray that should be traced from
// the camera through the specified pixel.
vec3 getRay(vec2 fragCoord, float windowWidth, float windowHeight, Camera camera){
    // Extract the x and y rotations of the camera
    float xRotation = camera.rotation.x;
    float yRotation = camera.rotation.y;
    
    // Create a matrix for rotating around the x-axis
    mat3 xRotationMatrix = mat3(
        1.0, 0.0, 0.0,
        0.0, cos(xRotation), sin(-xRotation),
        0.0, sin(xRotation), cos(xRotation)
    );
    
    // Create a matrix for rotating around the y-axis
    mat3 yRotationMatrix = mat3(
        cos(yRotation), 0.0, sin(yRotation),
        0.0, 1.0, 0.0,
        -sin(yRotation), 0.0, cos(yRotation)
    );
    
    // Calculate the position of the pixel on the projection plane
    vec3 pixelPosition = vec3(
        ((fragCoord.x/windowWidth)-0.5)*windowWidth/windowHeight,
        (fragCoord.y/windowHeight)-0.5,
        -camera.distanceFromProjectionPlane
    );
    
    // Rotate the pixel position according to the camera's orientation
    pixelPosition *= xRotationMatrix * yRotationMatrix;
    
    // Return the normalized direction of the ray from the camera through the pixel
    return normalize(pixelPosition);
}



// Calculates whether a given ray intersects a given sphere, and if it does, returns the
// time of the collision (i.e., the distance from the origin of the ray to the point of collision).
// If the ray does not intersect the sphere, returns -1.
float collideWith(Sphere sphere, vec3 rayOrigin, vec3 rayDirection){
    // Calculate b and c values for quadratic equation based on sphere's position,
    // ray's origin and direction, and sphere's radius
    float b = 2.0 * dot(rayDirection, rayOrigin - sphere.position);
    float c = pow(length(rayOrigin - sphere.position), 2.0) - pow(sphere.radius, 2.0);
    
    // Calculate discriminant (delta) of quadratic equation
    float delta = pow(b, 2.0) - 4.0 * c;
    
    // If delta is greater than 0, there is a solution (collision)
    if (delta > 0.0){
        // Calculate one solution to quadratic equation
        float t2 = (-b - sqrt(delta)) / 2.0;
        
        // Return the smaller of the two positive solutions (the other solution represents
        // a collision behind the ray's origin, which is not relevant in this case)
        if (t2 > 0.0){
            return t2;
        }
    }
    
    // If delta is less than or equal to 0, there is no solution (no collision)
    // or the smaller solution is negative (collision behind the ray's origin)
    return -1.0;
}



// Determines which sphere, if any, is the closest to the camera along a given ray.
// Outputs the index of the closest sphere and the point of collision.
// Returns true if a collision was found, false otherwise.
bool collideWithClosest(Sphere[SPHERE_COUNT] spheres, vec3 rayOrigin, vec3 rayDirection, out int oClosestSphereIndex, out vec3 oCollisionPoint){
    // Initialize variables to store index and collision time of closest sphere
    int closestSphereIndex;
    float timeOfClosestCollison = -1.0;
    
    // Iterate through all spheres
    for (int i = 0; i < SPHERE_COUNT; i++){
        // Check if ray collides with current sphere
        float time = collideWith(spheres[i], rayOrigin, rayDirection);
        
        // If there is a collision and it is the closest one so far,
        // update closestSphereIndex and timeOfClosestCollision
        if (time != -1.0 && (time < timeOfClosestCollison || timeOfClosestCollison == -1.0)){
            closestSphereIndex = i;
            timeOfClosestCollison = time;
        }
    }
    
    // If there was a collision with any sphere, set output variables and return true
    if (timeOfClosestCollison != -1.0){
        oClosestSphereIndex = closestSphereIndex;
        oCollisionPoint = rayOrigin + rayDirection*timeOfClosestCollison;
        return true;
    }
    
    // If there was no collision, return false
    return false;
}


vec3 getColor(Sphere[SPHERE_COUNT] spheres, vec3 rayOrigin, vec3 rayDirection, Camera camera, Light light){
    vec3 color = backgroundColor;
    vec3 _rayOrigin = rayOrigin;
    vec3 _rayDirection = rayDirection;
    float reflectiveness;
    
    for (int depth = 0; depth < DEPTH; depth++){
        int hitSphereIndex = -1;
        vec3 collisionPoint = vec3(0.0);
        bool rayCollided = collideWithClosest(spheres, _rayOrigin, _rayDirection, hitSphereIndex, collisionPoint);

        if (rayCollided){
            vec3 objectNormalVector = normalize(collisionPoint-spheres[hitSphereIndex].position);
            vec3 fromObjectToLightVector = normalize(light.position - collisionPoint);
            float diffuseStrength = max(0.0,dot(fromObjectToLightVector,objectNormalVector));

            vec3 fromObjectToCameraVector = normalize(camera.position - collisionPoint);
            float specular = pow(max(0.0,dot(objectNormalVector, normalize(fromObjectToLightVector+fromObjectToCameraVector))),spheres[hitSphereIndex].specularExponent);

            if (depth == 0){
                color = spheres[hitSphereIndex].color*light.color*(
                    spheres[hitSphereIndex].ambientStrength + diffuseStrength + spheres[hitSphereIndex].specularStrength*specular
                );
                reflectiveness = spheres[hitSphereIndex].reflectiveness;
            } else{
                color += reflectiveness*spheres[hitSphereIndex].color*light.color*(
                    spheres[hitSphereIndex].ambientStrength + diffuseStrength + spheres[hitSphereIndex].specularStrength*specular
                );
                reflectiveness*= spheres[hitSphereIndex].reflectiveness;
            }
            _rayOrigin = collisionPoint;
            _rayDirection = reflect(_rayDirection, objectNormalVector);
        }
    }
    
    return color;
}

vec3 render(vec2 fragCoord, float windowWidth, float windowHeight, Sphere[SPHERE_COUNT] spheres, Camera camera, Light light){
    vec3 color = vec3(0.0);
#ifdef ANTIALIASING
    for(float i = -1.0; i <= 1.0; i++) {
        for(float j = -1.0; j <= 1.0; j++) {
            vec3 rayDirection = getRay(
                fragCoord+vec2(i*0.3, j*0.3),
                windowWidth,
                windowHeight,
                camera
            );
    		color += getColor(spheres, camera.position, rayDirection, camera, light);
        }
    }
    color /= 9.0;
#else
    vec3 rayDirection = getRay(
        fragCoord,
        windowWidth,
        windowHeight,
        camera
    );
    color = getColor(spheres, camera.position, rayDirection, camera, light);
#endif
    return color;
}


float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float mouseInput = iMouse.x/100.0;
    float windowRatio = iResolution.y/iResolution.x;
    float fov = START_FOV+sin(iTime-mouseInput)*VARIANCE_FOV;
    
    Camera camera = Camera(
        vec3(sin(iTime-mouseInput)*5.0,CAMERA_HEIGHT+cos(iTime)/2.0,cos(iTime-mouseInput)*5.0),
        vec2(sin(iTime/2.0-mouseInput)/8.0-0.15,iTime-mouseInput),
        windowRatio,
        calculateDistanceFromProjectionPlane(windowRatio, fov)
    );
    
    float JUMP_HEIGHT = rand(vec2(1,3));
    
    float JUMP_RAND = rand(vec2(2,6));
    
    float BOUNCINESS_SMALL = abs(JUMP_HEIGHT * sin(iTime));
    float BOUNCINESS_BIG = abs(JUMP_HEIGHT * sin(iTime*1.0/JUMP_RAND));
    float BOUNCINESS_MEDIUM_SMALL = abs(JUMP_HEIGHT * sin(iTime*2.0));
    float BOUNCINESS_MEDIUM_BIG = abs(JUMP_HEIGHT * sin(iTime*1.5));
    
    float BIG_SPHERE_RADIUS = 1.0;
    float SMALL_SPHERE_RADIUS = 0.3;
    float MEDIUM_SMALL_SPHERE_RADIUS = 0.4;
    float MEDIUM_BIG_SPHERE_RADIUS = 0.7;
    
    
    Sphere spheres[SPHERE_COUNT] = Sphere[SPHERE_COUNT](
        Sphere( // Bigger Sphere
            vec3(-1.0, BOUNCINESS_BIG + BIG_SPHERE_RADIUS,0.0),    // position
            BIG_SPHERE_RADIUS,                    // radius 
            vec3(0.1,0.3,0.7),      // Color
            0.1,                    // Strength of the ambient lighting on the sphere
            0.7,                    // Strength of the specular lighting on the sphere
            128.0,                  // Exponent for the specular lighting on the sphere
            SPHERE_REFLECTIVENESS   // Reflectiveness of the sphere
        ),
        Sphere( // Small Sphere
            vec3(1.0,BOUNCINESS_SMALL + SMALL_SPHERE_RADIUS,0.0),
            SMALL_SPHERE_RADIUS,
            vec3(0.5,0.1,0.4),
            0.1,
            0.7,
            32.0,
            SPHERE_REFLECTIVENESS
        ),
        Sphere( // medium small Sphere
            vec3(-4.0,BOUNCINESS_MEDIUM_SMALL + MEDIUM_SMALL_SPHERE_RADIUS,2.0),
            MEDIUM_SMALL_SPHERE_RADIUS,
            vec3(0.2,0.8,0.2),
            0.2,
            0.7,
            32.0,
            SPHERE_REFLECTIVENESS
        ),
        Sphere( // medium big Sphere
            vec3(0.5,BOUNCINESS_MEDIUM_BIG + MEDIUM_BIG_SPHERE_RADIUS,4.3),
            MEDIUM_BIG_SPHERE_RADIUS,
            vec3(1.0,0.6,0.0),
            0.5,
            0.7,
            32.0,
            SPHERE_REFLECTIVENESS
        ),
        Sphere( // Floor
            vec3(0.0,-9000.0,0.0),
            9000.0,
            vec3(0.9,0.9,0.9),
            0.1,
            0.8,
            128.0,
            FLOOR_REFLECTIVENESS
        ),
        Sphere( // Ceiling
            vec3(0.0,9005.0,0.0),
            9000.0,
            vec3(0.9,0.9,0.9),
            0.1,
            0.8,
            128.0,
            CEILING_REFLECTIVENESS
        ),
        Sphere( // Left Wall 
            vec3(-9005.0,0.0,0.0),
            9000.0,
            vec3(1.0, 0.0, 0.0),
            0.1,
            0.8,
            128.0,
            WALL_REFLECTIVENESS
        ),
        Sphere( // Right Wall
            vec3(9005.0,0.0,0.0),
            9000.0,
            vec3(0.0, 0.0, 1.0),
            0.1,
            0.8,
            128.0,
            WALL_REFLECTIVENESS
        ),
        Sphere( // Front Wall
            vec3(0.0,0.0,-9005.0),
            9000.0,
            vec3(0.9,0.9,0.9),
            0.1,
            0.8,
            128.0,
            WALL_REFLECTIVENESS
        ),
        Sphere( // Back Wall
            vec3(0.0,0.0,+9015.0),
            9000.0,
            vec3(0.9,0.9,0.9),
            0.1,
            0.8,
            128.0,
            WALL_REFLECTIVENESS
        )
     );
     
     Light light = Light(
         vec3(3.0,3.0,6.0), // position?
         vec3(0.5,0.5,0.2),   // intensität?
         1.0
     );

    fragColor = vec4(render(fragCoord, iResolution.x, iResolution.y, spheres, camera, light),1.0);
}