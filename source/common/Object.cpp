//////////////////////////////////////////////////////////////////////////////
//
//  --- Object.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "common.h"
#include "Object.h"


using namespace std;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Sphere::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;
    result.t = raySphereIntersection(p0, V);
    if (result.t != std::numeric_limits<double>::infinity()) {
        result.P = p0 + (result.t * V);
        result.N = normalize(result.P - center);
    }
    return result;
}

/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin center and radius radius    : Find t ------- */
double Sphere::raySphereIntersection(const vec4 &p0, const vec4 &V) {
    vec2 t1t2 = vec2(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());

    //TODO: Ray-sphere intersection

    // t²d.d + 2td.(o-center) + ||o-center||²-r² = 0
    vec3 rayOrigin = vec3(p0.x, p0.y, p0.z); // o
    vec3 rayDirection = vec3(V.x, V.y, V.z); // d
    vec3 centerToRayOrigin = rayOrigin - center; // o-center
    // ||o-center||
    double normalizedCenterToRayOrigin = sqrt(
            pow(centerToRayOrigin.x, 2) + pow(centerToRayOrigin.y, 2) + pow(centerToRayOrigin.z, 2));

    double a = dot(rayDirection, rayDirection);
    double b = 2 * dot(rayDirection, centerToRayOrigin);
    double c = pow(normalizedCenterToRayOrigin, 2) - pow(radius, 2);
    double delta = pow(b, 2) - (4 * a * c);

    if (delta < 0) // No intersection, let t be infinity() and return it
        return std::numeric_limits<double>::infinity();

    if (delta == 0) // Only 1 intersection
        return (-b) / (2 * a);

    // Two intersections, keep the positive minimum between s1 and s2
    double s1 = (-b - sqrt(delta)) / (2 * a);
    double s2 = (-b + sqrt(delta)) / (2 * a);

    if (s1 < EPSILON && s2 < EPSILON)
        return std::numeric_limits<double>::infinity();
    else if (s1 < EPSILON)
        return s2;
    else if (s2 < EPSILON)
        return s1;
    else
        return min(s1, s2);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Square::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;
    result.t = raySquareIntersection(p0, V);
    if (result.t != std::numeric_limits<double>::infinity()) {
        result.P = p0 + result.t * V;
        result.N = vec4(normal, 0);
    }
    return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(const vec4 &p0, const vec4 &V) {
    vec3 p0_vec3 = vec3(p0.x, p0.y, p0.z);
    vec3 V_vec3 = vec3(V.x, V.y, V.z);
    vec3 point_vec3 = vec3(point.x, point.y, point.z);

    double denom = dot(V_vec3, normal);
    if (denom != 0) { // Avoid +infinity results when dividing by 0
        double t = (dot(point_vec3, normal) - dot(p0_vec3, normal)) / denom;
        if (t > EPSILON) { // The point is on a plan which is in front of the camera
            // Get the point coordinates on the plane
            vec3 P = p0_vec3 + t * V_vec3;

            // Get the mesh rectangle coordinates
            vec3 v0 = vec3(mesh.vertices[5].x, mesh.vertices[5].y, mesh.vertices[5].z); // Top left
            vec3 v1 = vec3(mesh.vertices[1].x, mesh.vertices[1].y, mesh.vertices[1].z); // Top right
            vec3 v2 = vec3(mesh.vertices[2].x, mesh.vertices[2].y, mesh.vertices[2].z); // Bottom right
            vec3 v3 = vec3(mesh.vertices[3].x, mesh.vertices[3].y, mesh.vertices[3].z); // Bottom left

            // Calculate the scalar products
            double value1 = dot(normal, cross(v1 - v0, P - v0));
            double value2 = dot(normal, cross(v2 - v1, P - v1));
            double value3 = dot(normal, cross(v3 - v2, P - v2));
            double value4 = dot(normal, cross(v0 - v3, P - v3));

            // If all the points are one the same side (all values are positive or negative), then the point is within the rectangle
            if (value1 >= 0 && value2 >= 0 && value3 >= 0 && value4 >= 0)
                return t;
            if (value1 <= 0 && value2 <= 0 && value3 <= 0 && value4 <= 0)
                return t;
        }
    }
    return std::numeric_limits<double>::infinity();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
bool Object::raySquareIntersection(const vec4 &p0, const vec4 &V, const vec3 &v0, const vec3 &v1, const vec3 &v2,
                                   const vec3 &v3) {
    vec3 p0_vec3 = vec3(p0.x, p0.y, p0.z);
    vec3 V_vec3 = vec3(V.x, V.y, V.z);

    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 normal = cross(edge1, edge2);
    if (dot(V, normal) > 0)
        normal = -normal;

    double denom = dot(V_vec3, normal);
    if (denom != 0) { // Avoid +infinity results when dividing by 0
        double t = (dot(v0, normal) - dot(p0_vec3, normal)) / denom;
        if (t > EPSILON) { // The point is on a plan which is in front of the camera
            // Get the point coordinates on the plane
            vec3 P = p0_vec3 + t * V_vec3;

            // Calculate the scalar products
            double value1 = dot(normal, cross(v1 - v0, P - v0));
            double value2 = dot(normal, cross(v2 - v1, P - v1));
            double value3 = dot(normal, cross(v3 - v2, P - v2));
            double value4 = dot(normal, cross(v0 - v3, P - v3));

            // If all the points are one the same side (all values are positive or negative), then the point is within the rectangle
            if (value1 >= 0 && value2 >= 0 && value3 >= 0 && value4 >= 0)
                return true;
            if (value1 <= 0 && value2 <= 0 && value3 <= 0 && value4 <= 0)
                return true;
        }
    }
    return false;
}

bool MeshObject::rayBoundingBoxIntersection(const vec4 &p0, const vec4 &V) {
    vec3 v0 = mesh.box_max; // Front top right
    vec3 v1 = vec3(mesh.box_max.x, mesh.box_min.y, mesh.box_max.z); // Front bottom right
    vec3 v2 = vec3(mesh.box_min.x, mesh.box_min.y, mesh.box_max.z); // Front bottom left
    vec3 v3 = vec3(mesh.box_min.x, mesh.box_max.y, mesh.box_max.z); // Front top left

    vec3 v4 = mesh.box_min; // Back bottom left
    vec3 v5 = vec3(mesh.box_min.x, mesh.box_max.y, mesh.box_min.z); // Back top left
    vec3 v6 = vec3(mesh.box_max.x, mesh.box_max.y, mesh.box_min.z); // Back top right
    vec3 v7 = vec3(mesh.box_max.x, mesh.box_min.y, mesh.box_min.z); // Back bottom right

    if (Object::raySquareIntersection(p0, V, v0, v1, v2, v3)) // Front face
        return true;
    if (Object::raySquareIntersection(p0, V, v4, v5, v6, v7)) // Back face
        return true;
    if (Object::raySquareIntersection(p0, V, v2, v3, v5, v4)) // Left face
        return true;
    if (Object::raySquareIntersection(p0, V, v0, v1, v7, v6)) // Right face
        return true;
    if (Object::raySquareIntersection(p0, V, v0, v3, v5, v6)) // Top face
        return true;
    return Object::raySquareIntersection(p0, V, v1, v2, v4, v7); // Bottom face
}

Object::IntersectionValues MeshObject::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;
    result.t = numeric_limits<double>::infinity();
    if (rayBoundingBoxIntersection(p0, V))
        rayMeshObjectIntersection(p0, V, &result);
    return result;
}

void MeshObject::rayTriangleIntersection(const vec3 &origin, const vec3 &ray, const vec3 &v0, const vec3 &v1,
                                         const vec3 &v2, const vec3 &n0, const vec3 &n1, const vec3 &n2,
                                         IntersectionValues *result) {
    const float epsilon = 0.0000001f;
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;

    vec3 h = cross(ray, edge2);
    float a = dot(edge1, h);
    if (a > -epsilon && a < epsilon) // Ray parallel to triangle
        return;

    float f = 1.0f / a;
    vec3 s = origin - v0;
    float u = f * dot(s, h);
    if (u < 0 || u > 1)
        return;

    vec3 q = cross(s, edge1);
    float v = f * dot(ray, q);
    if (v < 0 || u + v > 1)
        return;

    // There is an intersection, calculate if it is in front of the camera
    float t = f * dot(edge2, q);
    if (t <= epsilon) // Behind
        return;

    // Intersection in front of the camera, change the results if it is the closest to the camera
    if (t < result->t) {
        result->t = t;
        result->P = vec4(origin + (result->t * ray), 1);

        vec3 N = normalize(cross(edge1, edge2));
        if (dot(N, ray) > 0)
            result->N = vec4(-N, 0);
        else
            result->N = vec4(N, 0);
    }
}

void MeshObject::rayMeshObjectIntersection(const vec4 &p0, const vec4 &V, IntersectionValues *result) {
    vec3 origin = vec3(p0.x, p0.y, p0.z);
    vec3 ray = vec3(V.x, V.y, V.z);

    for (int index = 0; index < mesh.getNumTri(); index++) {
        int i = index * 3;

        vec3 v0 = vec3(mesh.vertices[i].x, mesh.vertices[i].y, mesh.vertices[i].z);
        vec3 n0 = vec3(mesh.normals[i].x, mesh.normals[i].y, mesh.normals[i].z);

        vec3 v1 = vec3(mesh.vertices[i + 1].x, mesh.vertices[i + 1].y, mesh.vertices[i + 1].z);
        vec3 n1 = vec3(mesh.normals[i + 1].x, mesh.normals[i + 1].y, mesh.normals[i + 1].z);

        vec3 v2 = vec3(mesh.vertices[i + 2].x, mesh.vertices[i + 2].y, mesh.vertices[i + 2].z);
        vec3 n2 = vec3(mesh.normals[i + 2].x, mesh.normals[i + 2].y, mesh.normals[i + 2].z);

        rayTriangleIntersection(origin, ray, v0, v1, v2, n0, n1, n2, result);
    }
}
