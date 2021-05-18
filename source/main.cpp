//////////////////////////////////////////////////////////////////////////////
//
//  --- main.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "SourcePath.h"
#include <io.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <Chrono.h>

using namespace Angel;
using namespace std;
using namespace std::chrono;

typedef vec4 color4;
typedef vec4 point4;


//Scene variables
enum {
    _SPHERE, _SQUARE, _BOX
};
int scene = _SPHERE; //Simple sphere, square or cornell box
std::vector<Object *> sceneObjects;
point4 lightPosition;
color4 lightColor;
point4 cameraPosition;

vector<point4> lightPoints;
double lightSize = 0.5;
int nbPoints = 50;

//Recursion depth for raytracer
int maxDepth = 5;

void initGL();
vec4 castRay(const vec4& p0, const vec4& dir, Object *lastHitObject, int depth, float currentKr, bool debug);

namespace GLState {
    int window_width, window_height;

    bool render_line;

    std::vector<GLuint> objectVao;
    std::vector<GLuint> objectBuffer;

    GLuint vPosition, vNormal, vTexCoord;

    GLuint program;

    // Model-view and projection matrices uniform location
    GLuint ModelView, ModelViewLight, NormalMatrix, Projection;

    //==========Trackball Variables==========
    static float curquat[4], lastquat[4];
    /* current transformation matrix */
    static float curmat[4][4];
    mat4 curmat_a;
    /* actual operation  */
    static int scaling;
    static int moving;
    static int panning;
    /* starting "moving" coordinates */
    static int beginx, beginy;
    /* ortho */
    float ortho_x, ortho_y;
    /* current scale factor */
    static float scalefactor;

    mat4 projection;
    mat4 sceneModelView;

    color4 light_ambient;
    color4 light_diffuse;
    color4 light_specular;

};

/* ------------------------------------------------------- */
/* -- PNG receptor class for use with pngdecode library -- */
class rayTraceReceptor : public cmps3120::png_receptor {
private:
    const unsigned char *buffer;
    unsigned int width;
    unsigned int height;
    int channels;

public:
    rayTraceReceptor(const unsigned char *use_buffer,
                     unsigned int width,
                     unsigned int height,
                     int channels) {
        this->buffer = use_buffer;
        this->width = width;
        this->height = height;
        this->channels = channels;
    }

    cmps3120::png_header get_header() {
        cmps3120::png_header header;
        header.width = width;
        header.height = height;
        header.bit_depth = 8;
        switch (channels) {
            case 1:
                header.color_type = cmps3120::PNG_GRAYSCALE;
                break;
            case 2:
                header.color_type = cmps3120::PNG_GRAYSCALE_ALPHA;
                break;
            case 3:
                header.color_type = cmps3120::PNG_RGB;
                break;
            default:
                header.color_type = cmps3120::PNG_RGBA;
                break;
        }
        return header;
    }

    cmps3120::png_pixel get_pixel(unsigned int x, unsigned int y, unsigned int level) {
        cmps3120::png_pixel pixel;
        unsigned int idx = y * width + x;
        /* pngdecode wants 16-bit color values */
        pixel.r = buffer[4 * idx] * 257;
        pixel.g = buffer[4 * idx + 1] * 257;
        pixel.b = buffer[4 * idx + 2] * 257;
        pixel.a = buffer[4 * idx + 3] * 257;
        return pixel;
    }
};

/* -------------------------------------------------------------------------- */
/* ----------------------  Write Image to Disk  ----------------------------- */
bool write_image(const char *filename, const unsigned char *Src,
                 int Width, int Height, int channels) {
    cmps3120::png_encoder the_encoder;
    cmps3120::png_error result;
    rayTraceReceptor image(Src, Width, Height, channels);
    the_encoder.set_receptor(&image);
    result = the_encoder.write_file(filename);
    if (result == cmps3120::PNG_DONE)
        std::cerr << "finished writing " << filename << "." << std::endl;
    else
        std::cerr << "write to " << filename << " returned error code " << result << "." << std::endl;
    return result == cmps3120::PNG_DONE;
}


/* -------------------------------------------------------------------------- */
/* -------- Given OpenGL matrices find ray in world coordinates of ---------- */
/* -------- window position x,y --------------------------------------------- */
std::vector<vec4> findRay(GLdouble x, GLdouble y) {

    y = GLState::window_height - y;

    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    GLdouble modelViewMatrix[16];
    GLdouble projectionMatrix[16];
    for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = 0; j < 4; j++) {
            modelViewMatrix[j * 4 + i] = GLState::sceneModelView[i][j];
            projectionMatrix[j * 4 + i] = GLState::projection[i][j];
        }
    }


    GLdouble nearPlaneLocation[3];
    _gluUnProject(x, y, 0.0, modelViewMatrix, projectionMatrix,
                  viewport, &nearPlaneLocation[0], &nearPlaneLocation[1],
                  &nearPlaneLocation[2]);

    GLdouble farPlaneLocation[3];
    _gluUnProject(x, y, 1.0, modelViewMatrix, projectionMatrix,
                  viewport, &farPlaneLocation[0], &farPlaneLocation[1],
                  &farPlaneLocation[2]);


    vec4 ray_origin = vec4(nearPlaneLocation[0], nearPlaneLocation[1], nearPlaneLocation[2], 1.0);
    vec3 temp = vec3(farPlaneLocation[0] - nearPlaneLocation[0],
                     farPlaneLocation[1] - nearPlaneLocation[1],
                     farPlaneLocation[2] - nearPlaneLocation[2]);
    temp = normalize(temp);
    vec4 ray_dir = vec4(temp.x, temp.y, temp.z, 0.0);

    std::vector<vec4> result(2);
    result[0] = ray_origin;
    result[1] = ray_dir;

    return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
bool intersectionSort(Object::IntersectionValues i, Object::IntersectionValues j) {
    return (i.t < j.t);
}

color4 calculateAmbientColor(const Object::IntersectionValues& intersectionValue) {
    Object::ShadingValues shVal = sceneObjects[intersectionValue.ID_]->shadingValues;

    // Ambient light intensity
    color4 ambientMaterial = shVal.color * shVal.Ka; ambientMaterial.w = 1;
    color4 ambient = GLState::light_ambient * ambientMaterial; ambient.w = 1;

    return ambient;
}

color4 calculateIllumination(const Object::IntersectionValues& intersectionValue, const vec4& dir) {

    vec4 P = intersectionValue.P; P.w = 1;
    vec4 N = normalize(intersectionValue.N); N.w = 0;
    vec4 L = normalize(lightPosition - P); L.w = 0;
    vec4 V = normalize(-dir); V.w = 0;
    vec4 R = normalize(-reflect(L, N));

    Object::ShadingValues shVal = sceneObjects[intersectionValue.ID_]->shadingValues;

    // Ambient light intensity
    color4 ambientMaterial = shVal.color * shVal.Ka; ambientMaterial.w = 1;
    color4 ambient = GLState::light_ambient * ambientMaterial; ambient.w = 1;

    // Diffuse light intensity
    color4 diffuseMaterial = shVal.color * shVal.Kd; diffuseMaterial.w = 1;
    float diffuseFactor = max(dot(L, N), 0.0f);
    color4 diffuse = GLState::light_diffuse * diffuseMaterial * diffuseFactor; diffuse.w = 1;

    // Specular light intensity
    vec4 specularMaterial = shVal.Ks;
    float specularFactor = pow(max(dot(V, R), 0.0f), shVal.Kn);
    color4 specular = GLState::light_specular * specularMaterial * specularFactor; specular.w = 1; //  * shVal.color

    color4 final = ambient + diffuse + specular;

    final.limit(); // Best of the 3 methods

    final.w = 1;
    return final;
}

// Cast a ray to the light source in the scene, if the ray intersects an object and the intersection occurs before
// the ray arrived to the light source, then the light from the light source is blocked by an object
// If the light is blocked return true, else return false
float shadowFactorTransparency(const vec4 &point, const vec4 &lightPos = lightPosition) {
    vec4 pointToLight = lightPos - point; pointToLight.w = 0;

    vec4 pointToLightNormalized = normalize(pointToLight); pointToLightNormalized.w = 0;
    vec4 epsilonPoint = point + EPSILON * pointToLightNormalized; epsilonPoint.w = 1;

    vec4 epsilonPointToLight = lightPos - epsilonPoint; epsilonPointToLight.w = 0;
    vec4 epsilonPointToLightNormalized = normalize(epsilonPointToLight); epsilonPointToLightNormalized.w = 0;

    std::vector<Object::IntersectionValues> intersectionValuesVector;
    for (int i = 0; i < sceneObjects.size(); i++) {
        intersectionValuesVector.push_back(sceneObjects[i]->intersect(epsilonPoint, epsilonPointToLightNormalized));
        intersectionValuesVector[intersectionValuesVector.size() - 1].ID_ = i;
    }

    float coef = 1;
    for (auto &intersectionValues : intersectionValuesVector) {
        double lengthFromEpsilonPointToLight = length(epsilonPointToLight);

        if (intersectionValues.t != std::numeric_limits<double>::infinity() && intersectionValues.t >= 0) {
            vec4 epsilonPointToP = intersectionValues.P - epsilonPoint;
            double lengthFromEpsilonPointToP = length(epsilonPointToP);
            if (lengthFromEpsilonPointToP > 0 && lengthFromEpsilonPointToP < lengthFromEpsilonPointToLight) { // Point is in shadow
                Object::ShadingValues shVal = sceneObjects[intersectionValues.ID_]->shadingValues;
                coef *= 0.9f * shVal.Kt;
            }
        }
    }

    return coef;
}

double randomDouble(double min, double max) {
    double d = (double) rand() / RAND_MAX;
    return min + d * (max - min);
}

void generateLightPoints() {
    lightPoints.clear();

    uint64_t ns = duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
    srand(ns);

    int steps = 2;
    double maxStepSize = lightSize / (2 * steps);

    for (int i = 0; i < nbPoints; i++) {
        double x = lightPosition.x;
        double z = lightPosition.z;
        for (int j = 0; j < steps; j++) {
            x += randomDouble(-maxStepSize, +maxStepSize);
            z += randomDouble(-maxStepSize, +maxStepSize);
        }
        lightPoints.emplace_back(x, lightPosition.y, z, lightPosition.w);
    }
}

float factorPointIsInShadow(const vec4 &point) {
    float total = 0;
    for (point4 &lightPoint : lightPoints)
        total += shadowFactorTransparency(point, lightPoint);
    return total / lightPoints.size();
}

static color4 maxColor(const color4& color1, const color4& color2) {
    if (color1.x > color2.x) return color1;
    if (color1.x < color2.x) return color2;
    if (color1.y > color2.y) return color1;
    if (color1.y < color2.y) return color2;
    if (color1.z > color2.z) return color1;
    if (color1.z < color2.z) return color2;
    if (color1.w > color2.w) return color1;
    if (color1.w < color2.w) return color2;
    return color1; // Same colors
}

static color4 interpolate(const color4& color1, const color4& color2, double factor) {
    if (maxColor(color1, color2) == color1)
        return color1;

    return (1 - factor) * color1 + factor * color2;
}

color4 getRefractedColor(const vec4& dir, const vec4& normal, const point4& newP0, const Object::ShadingValues& shVal,
                         color4 objectColor, Object *lastHitObject, float currentKr, int depth, bool exitingObject, bool debug) {
    objectColor *= (1 - shVal.Kt);

    vec4 refractedNormal = normal;
    float newKr = shVal.Kr;
    if (exitingObject) {
        newKr = 1.0f;
        refractedNormal = -refractedNormal;
    }

    vec4 refractedRay = refract(dir, refractedNormal, currentKr, newKr);
    if (debug) {
        cout << "Transparency = " << shVal.Kt << ", current coef. = " << currentKr << ", refraction coef. = " << newKr << endl;
        cout << "Incident ray = " << dir << endl;
        cout << "Surface normal = " << refractedNormal << endl;
        cout << "Refracted ray = " << refractedRay << endl;
    }

    color4 refractedObjectColor = shVal.Kt * castRay(newP0, refractedRay, lastHitObject, depth + 1, newKr, debug);

    return objectColor + refractedObjectColor;
}

/* -------------------------------------------------------------------------- */
/* ----------  cast Ray = p0 + t*dir and intersect with sphere      --------- */
/* ----------  return color, right now shading is approx based      --------- */
/* ----------  depth                                                --------- */
vec4 castRay(const vec4& p0, const vec4& dir, Object *lastHitObject, int depth, float currentKr, bool debug) {
    color4 color = vec4(0.0, 0.0, 0.0, 0.0);

    std::vector<Object::IntersectionValues> intersectionValuesVector;
    for (int i = 0; i < sceneObjects.size(); i++) {
        intersectionValuesVector.push_back(sceneObjects[i]->intersect(p0, dir));
        intersectionValuesVector[intersectionValuesVector.size() - 1].ID_ = i;
    }

    double min = std::numeric_limits<double>::infinity();
    int id = -1;
    for (auto &intersectionValues : intersectionValuesVector) {
        if (intersectionValues.t != std::numeric_limits<double>::infinity() && intersectionValues.t < min) {
            if (depth == 0 or intersectionValues.t > EPSILON) {
                min = intersectionValues.t;
                id = intersectionValues.ID_;
            }
        }
    }

    if (id == -1) { // No object hit
        return color;
    }

    // Lighting and shadows
    bool exitingObject = false;
    if (lastHitObject == sceneObjects[id])
        exitingObject = true;

    lastHitObject = sceneObjects[id];
    Object::ShadingValues shVal = lastHitObject->shadingValues;
    Object::IntersectionValues iVal = intersectionValuesVector[id];
    if (debug) {
        cout << "-----" << endl;
        cout << "Hit: " << lastHitObject->name << " at " << iVal.P << endl;
    }

    color4 ambientColor = calculateAmbientColor(intersectionValuesVector[id]);
    color4 phongColor = calculateIllumination(intersectionValuesVector[id], dir);

    double shadowFactor = factorPointIsInShadow(intersectionValuesVector[id].P);
    color4 objectColor = interpolate(ambientColor, phongColor, shadowFactor);

    if (depth == maxDepth) {
        color = objectColor;
    }
    else {
        point4 newP0 = iVal.P;
        vec4 normal = iVal.N;

        if (shVal.Ks == 0) { // Non-mirror object
            if (shVal.Kt == 0) { // Non-transparent AND non-transparent object
                color = objectColor;
            }
            else { // Non-mirror but transparent object
                color = getRefractedColor(dir, normal, newP0, shVal, objectColor, lastHitObject, currentKr, depth, exitingObject, debug);
            }
        }
        else { // Mirror object
            if (shVal.Kt != 0) { // Mirror AND transparent object
                objectColor = getRefractedColor(dir, normal, newP0, shVal, objectColor, lastHitObject, currentKr, depth, exitingObject, debug);
            }
//            color4 objectColor = interpolate(ambientColor, phongColor, shadowFactor);
            objectColor *= (1 - shVal.Ks);

            vec4 reflectedRay = normalize(reflect(dir, iVal.N));

            color4 reflectedObjetColor = shVal.Ks * castRay(newP0, reflectedRay, lastHitObject, depth + 1, currentKr, debug);

            color = objectColor + reflectedObjetColor;
        }
    }

    color.w = 1.0;
//    cout << "Color = " << color << endl;
    return color;
}

/* -------------------------------------------------------------------------- */
/* ---------  Some debugging code: cast Ray = p0 + t*dir  ------------------- */
/* ---------  and print out what it hits =                ------------------- */
void castRayDebug(const vec4& p0, const vec4& dir) {

    cout << "### DEBUG START ###" << endl;
    castRay(p0, dir, nullptr, 0, 1.0f, true);
    cout << "# DEBUG END #" << endl << endl << endl << endl << endl << endl;

//    std::vector<Object::IntersectionValues> intersections;
//
//    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
//        intersections.push_back(sceneObjects[i]->intersect(p0, dir));
//        intersections[intersections.size() - 1].ID_ = i;
//    }
//
//    for (unsigned int i = 0; i < intersections.size(); i++) {
//        if (intersections[i].t != std::numeric_limits<double>::infinity()) {
//            std::cout << "Hit " << intersections[i].name << " " << intersections[i].ID_ << endl;
//            std::cout << "P: " << intersections[i].P << endl;
//            std::cout << "N: " << intersections[i].N << endl;
//            vec4 L = lightPosition - intersections[i].P;
//            L = normalize(L);
//            std::cout << "L: " << L << endl;
//        }
//    }

}

/* -------------------------------------------------------------------------- */
/* ------------  Ray trace our scene.  Output color to image and    --------- */
/* -----------   Output color to image and save to disk             --------- */
void rayTrace() {
    auto* chrono = new Chrono("Render time");
    cout << "Starting the rendering..." << endl;

    auto *buffer = new unsigned char[GLState::window_width * GLState::window_height * 4];
    generateLightPoints();

    auto* findRayChrono = new Chrono("findRay total time");
    findRayChrono->pause();
    auto* castRayChrono = new Chrono("castRay total time");
    castRayChrono->pause();

    double currentPx = 0;
    double totalPx = GLState::window_width * GLState::window_height;
    double currentProgressionPercentage = 0;
    double printEveryPercentage = 5;

    for (unsigned int i = 0; i < GLState::window_width; i++) {
        for (unsigned int j = 0; j < GLState::window_height; j++) {
            unsigned int idx = j * GLState::window_width + i;
            findRayChrono->unpause();
            std::vector<vec4> ray_o_dir = findRay(i, j);
            findRayChrono->pause();

            castRayChrono->unpause();
            vec4 color = castRay(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z, 0.0), nullptr, 0, 1.0f, false);
            castRayChrono->pause();
            buffer[4 * idx] = color.x * 255;
            buffer[4 * idx + 1] = color.y * 255;
            buffer[4 * idx + 2] = color.z * 255;
            buffer[4 * idx + 3] = color.w * 255;

            currentPx++;
            if ((currentPx / totalPx) * 100 >= currentProgressionPercentage + printEveryPercentage) {
                currentProgressionPercentage += printEveryPercentage;
                cout << round(currentProgressionPercentage) << "% (ETR = " << round(
                        chrono->getETR(currentProgressionPercentage)) << "s)" << endl;
            }
        }
    }

    findRayChrono->displayCurrentDuration();
    castRayChrono->displayCurrentDuration();

    auto* writeChrono = new Chrono("Write image");
    write_image("output.png", buffer, GLState::window_width, GLState::window_height, 4);
    delete[] buffer;
    writeChrono->displayCurrentDuration();

    cout << "Rendering ended." << endl;
    chrono->displayCurrentDuration();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void error_callback(int error, const char *description) {
    fprintf(stderr, "Error: %s\n", description);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initCornellBox() {
    cameraPosition = point4(0.0, 0.0, 6.0, 1.0);
    lightPosition = point4(0.0, 1.5, 0.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    float ka = 0.1;
    float kd = 0.9;
    float ks = 0.0;
    float kn = 16.0;
    float kt = 0.0;
    float kr = 1.0; // Refraction coef. ALWAYS >= 1

    { //Back Wall
        sceneObjects.push_back(new Square("Back Wall", Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.5, 1.0, 1.0, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.9f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    { //Left Wall
        sceneObjects.push_back(new Square("Left Wall", RotateY(90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.8, 0.0, 0.0, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.0f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    { //Right Wall
        sceneObjects.push_back(
                new Square("Right Wall", RotateY(-90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.5, 0.0, 0.5, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.0f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    { //Floor
        sceneObjects.push_back(new Square("Floor", RotateX(-90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.3, 1.0, 0.3, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.0f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    { //Ceiling
        sceneObjects.push_back(new Square("Ceiling", RotateX(90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.5, 0.0, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.5f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    { //Front Wall
        sceneObjects.push_back(
                new Square("Front Wall", RotateY(180) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.2, 0.2, 1.0, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.0f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }


    {
        sceneObjects.push_back(new Sphere("Mirrored sphere 1", vec3(-1.0, -0.5, -1.0), 0.75));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.1, 0.5, 0.1, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 1.0f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    {
        sceneObjects.push_back(new Sphere("Mirror sphere 2", vec3(1, 1.25, -1), 0.5));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.9, 0.6, 0.1, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.75f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.0f;
        _shadingValues.Kr = kr + 0.0f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    {
        sceneObjects.push_back(new Sphere("Glass sphere", vec3(1.0, -1.25, 0.5), 0.75));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.5, 0.5, 1.0);
        _shadingValues.Ka = ka + 0.0f;
        _shadingValues.Kd = kd + 0.0f;
        _shadingValues.Ks = ks + 0.1f;
        _shadingValues.Kn = kn + 0.0f;
        _shadingValues.Kt = kt + 0.9f;
        _shadingValues.Kr = kr + 0.5f;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSphere() {
    cameraPosition = point4(0.0, 0.0, 3.0, 1.0);
    lightPosition = point4(0.0, 0.0, 4.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    {
        sceneObjects.push_back(new Sphere("Diffuse sphere", vec3(1.0, 0.0, 0.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0, 0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    {
        sceneObjects.push_back(new Sphere("Diffuse sphere 2", vec3(-1.0, 0.0, 0.0)));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0, 1, 0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSquare() {
    cameraPosition = point4(0.0, 0.0, 3.0, 1.0);
    lightPosition = point4(0.0, 0.0, 4.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    { //Back Wall
        sceneObjects.push_back(new Square("Unit Square"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0, 0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
        if (scene != _SPHERE) {
            initUnitSphere();
            initGL();
            scene = _SPHERE;
        }
    }
    if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
        if (scene != _SQUARE) {
            initUnitSquare();
            initGL();
            scene = _SQUARE;
        }
    }
    if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
        if (scene != _BOX) {
            initCornellBox();
            initGL();
            scene = _BOX;
        }
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS) {
        rayTrace();
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void mouseClick(GLFWwindow *window, int button, int action, int mods) {

    if (GLFW_RELEASE == action) {
        GLState::moving = GLState::scaling = GLState::panning = false;
        return;
    }

    if (mods & GLFW_MOD_SHIFT) {
        GLState::scaling = true;
    } else if (mods & GLFW_MOD_ALT) {
        GLState::panning = true;
    } else {
        GLState::moving = true;
        TrackBall::trackball(GLState::lastquat, 0, 0, 0, 0);
    }

    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    GLState::beginx = xpos;
    GLState::beginy = ypos;

    std::vector<vec4> ray_o_dir = findRay(xpos, ypos);
    castRayDebug(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z, 0.0));

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void mouseMove(GLFWwindow *window, double x, double y) {

    int W, H;
    glfwGetFramebufferSize(window, &W, &H);


    float dx = (x - GLState::beginx) / (float) W;
    float dy = (GLState::beginy - y) / (float) H;

    if (GLState::panning) {
        GLState::ortho_x += dx;
        GLState::ortho_y += dy;

        GLState::beginx = x;
        GLState::beginy = y;
        return;
    } else if (GLState::scaling) {
        GLState::scalefactor *= (1.0f + dx);

        GLState::beginx = x;
        GLState::beginy = y;
        return;
    } else if (GLState::moving) {
        TrackBall::trackball(GLState::lastquat,
                             (2.0f * GLState::beginx - W) / W,
                             (H - 2.0f * GLState::beginy) / H,
                             (2.0f * x - W) / W,
                             (H - 2.0f * y) / H
        );

        TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
        TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

        GLState::beginx = x;
        GLState::beginy = y;
        return;
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initGL() {

    GLState::light_ambient = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);
    GLState::light_diffuse = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);
    GLState::light_specular = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);


    std::string vshader = source_path + "/shaders/vshader.glsl";
    std::string fshader = source_path + "/shaders/fshader.glsl";

    GLchar *vertex_shader_source = readShaderSource(vshader.c_str());
    GLchar *fragment_shader_source = readShaderSource(fshader.c_str());

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, (const GLchar **) &vertex_shader_source, NULL);
    glCompileShader(vertex_shader);
    check_shader_compilation(vshader, vertex_shader);

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, (const GLchar **) &fragment_shader_source, NULL);
    glCompileShader(fragment_shader);
    check_shader_compilation(fshader, fragment_shader);

    GLState::program = glCreateProgram();
    glAttachShader(GLState::program, vertex_shader);
    glAttachShader(GLState::program, fragment_shader);

    glLinkProgram(GLState::program);
    check_program_link(GLState::program);

    glUseProgram(GLState::program);

    glBindFragDataLocation(GLState::program, 0, "fragColor");

    // set up vertex arrays
    GLState::vPosition = glGetAttribLocation(GLState::program, "vPosition");
    GLState::vNormal = glGetAttribLocation(GLState::program, "vNormal");

    // Retrieve transformation uniform variable locations
    GLState::ModelView = glGetUniformLocation(GLState::program, "ModelView");
    GLState::NormalMatrix = glGetUniformLocation(GLState::program, "NormalMatrix");
    GLState::ModelViewLight = glGetUniformLocation(GLState::program, "ModelViewLight");
    GLState::Projection = glGetUniformLocation(GLState::program, "Projection");

    GLState::objectVao.resize(sceneObjects.size());
    glGenVertexArrays(sceneObjects.size(), &GLState::objectVao[0]);

    GLState::objectBuffer.resize(sceneObjects.size());
    glGenBuffers(sceneObjects.size(), &GLState::objectBuffer[0]);

    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
        glBindVertexArray(GLState::objectVao[i]);
        glBindBuffer(GL_ARRAY_BUFFER, GLState::objectBuffer[i]);
        size_t vertices_bytes = sceneObjects[i]->mesh.vertices.size() * sizeof(vec4);
        size_t normals_bytes = sceneObjects[i]->mesh.normals.size() * sizeof(vec3);

        glBufferData(GL_ARRAY_BUFFER, vertices_bytes + normals_bytes, NULL, GL_STATIC_DRAW);
        size_t offset = 0;
        glBufferSubData(GL_ARRAY_BUFFER, offset, vertices_bytes, &sceneObjects[i]->mesh.vertices[0]);
        offset += vertices_bytes;
        glBufferSubData(GL_ARRAY_BUFFER, offset, normals_bytes, &sceneObjects[i]->mesh.normals[0]);

        glEnableVertexAttribArray(GLState::vNormal);
        glEnableVertexAttribArray(GLState::vPosition);

        glVertexAttribPointer(GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
        glVertexAttribPointer(GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(vertices_bytes));

    }


    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    glClearColor(0.8, 0.8, 1.0, 1.0);

    //Quaternion trackball variables, you can ignore
    GLState::scaling = 0;
    GLState::moving = 0;
    GLState::panning = 0;
    GLState::beginx = 0;
    GLState::beginy = 0;

    TrackBall::matident(GLState::curmat);
    TrackBall::trackball(GLState::curquat, 0.0f, 0.0f, 0.0f, 0.0f);
    TrackBall::trackball(GLState::lastquat, 0.0f, 0.0f, 0.0f, 0.0f);
    TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
    TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

    GLState::scalefactor = 1.0;
    GLState::render_line = false;

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void drawObject(Object *object, GLuint vao, GLuint buffer) {

    color4 material_ambient(object->shadingValues.color.x * object->shadingValues.Ka,
                            object->shadingValues.color.y * object->shadingValues.Ka,
                            object->shadingValues.color.z * object->shadingValues.Ka, 1.0);
    color4 material_diffuse(object->shadingValues.color.x,
                            object->shadingValues.color.y,
                            object->shadingValues.color.z, 1.0);
    color4 material_specular(object->shadingValues.Ks,
                             object->shadingValues.Ks,
                             object->shadingValues.Ks, 1.0);
    float material_shininess = object->shadingValues.Kn;

    color4 ambient_product = GLState::light_ambient * material_ambient;
    color4 diffuse_product = GLState::light_diffuse * material_diffuse;
    color4 specular_product = GLState::light_specular * material_specular;

    glUniform4fv(glGetUniformLocation(GLState::program, "AmbientProduct"), 1, ambient_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "DiffuseProduct"), 1, diffuse_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "SpecularProduct"), 1, specular_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "LightPosition"), 1, lightPosition);
    glUniform1f(glGetUniformLocation(GLState::program, "Shininess"), material_shininess);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glVertexAttribPointer(GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
    glVertexAttribPointer(GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0,
                          BUFFER_OFFSET(object->mesh.vertices.size() * sizeof(vec4)));

    mat4 objectModelView = GLState::sceneModelView * object->getModelView();


    glUniformMatrix4fv(GLState::ModelViewLight, 1, GL_TRUE, GLState::sceneModelView);
    glUniformMatrix3fv(GLState::NormalMatrix, 1, GL_TRUE, Normal(objectModelView));
    glUniformMatrix4fv(GLState::ModelView, 1, GL_TRUE, objectModelView);

    glDrawArrays(GL_TRIANGLES, 0, object->mesh.vertices.size());
}


int main(void) {

    GLFWwindow *window;

    glfwSetErrorCallback(error_callback);

    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_SAMPLES, 4);

    window = glfwCreateWindow(768, 768, "Raytracer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseClick);
    glfwSetCursorPosCallback(window, mouseMove);


    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
    glfwSwapInterval(1);

    switch (scene) {
        case _SPHERE:
            initUnitSphere();
            break;
        case _SQUARE:
            initUnitSquare();
            break;
        case _BOX:
            initCornellBox();
            break;
    }

    initGL();

    while (!glfwWindowShouldClose(window)) {

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);

        GLState::window_height = height;
        GLState::window_width = width;

        glViewport(0, 0, width, height);


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mat4 track_ball = mat4(GLState::curmat[0][0], GLState::curmat[1][0],
                               GLState::curmat[2][0], GLState::curmat[3][0],
                               GLState::curmat[0][1], GLState::curmat[1][1],
                               GLState::curmat[2][1], GLState::curmat[3][1],
                               GLState::curmat[0][2], GLState::curmat[1][2],
                               GLState::curmat[2][2], GLState::curmat[3][2],
                               GLState::curmat[0][3], GLState::curmat[1][3],
                               GLState::curmat[2][3], GLState::curmat[3][3]);

        GLState::sceneModelView = Translate(-cameraPosition) *   //Move Camera Back
                                  Translate(GLState::ortho_x, GLState::ortho_y, 0.0) *
                                  track_ball *                   //Rotate Camera
                                  Scale(GLState::scalefactor,
                                        GLState::scalefactor,
                                        GLState::scalefactor);   //User Scale

        GLfloat aspect = GLfloat(width) / height;

        switch (scene) {
            case _SPHERE:
            case _SQUARE:
                GLState::projection = Perspective(45.0, aspect, 0.01, 100.0);
                break;
            case _BOX:
                GLState::projection = Perspective(45.0, aspect, 4.5, 100.0);
                break;
        }

        glUniformMatrix4fv(GLState::Projection, 1, GL_TRUE, GLState::projection);

        for (unsigned int i = 0; i < sceneObjects.size(); i++) {
            drawObject(sceneObjects[i], GLState::objectVao[i], GLState::objectBuffer[i]);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
