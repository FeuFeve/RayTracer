//////////////////////////////////////////////////////////////////////////////
//
//  --- Object.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "common.h"

using namespace std;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Sphere::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;
    vec2 t1t2 = raySphereIntersection(p0, V);
    result.t = t1t2.x;
    result.t2 = t1t2.y;
    if (result.t != std::numeric_limits<double>::infinity()) {
        result.P = p0 + (result.t * V);
        result.P2 = p0 + (t1t2.y * V);
        result.N = normalize(result.P - center);
    }
    return result;
}

/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin center and radius radius    : Find t ------- */
vec2 Sphere::raySphereIntersection(const vec4& p0, const vec4& V) {
    vec2 t1t2 = vec2(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());

    //TODO: Ray-sphere intersection

    // t²d.d + 2td.(o-center) + ||o-center||²-r² = 0
    vec3 rayOrigin = vec3(p0.x, p0.y, p0.z); // o
    vec3 rayDirection = vec3(V.x, V.y, V.z); // d
    vec3 centerToRayOrigin = rayOrigin - center; // o-center
    // ||o-center||
    double normalizedCenterToRayOrigin = sqrt(pow(centerToRayOrigin.x, 2) + pow(centerToRayOrigin.y, 2) + pow(centerToRayOrigin.z, 2));

    double a = dot(rayDirection, rayDirection);
    double b = 2 * dot(rayDirection, centerToRayOrigin);
    double c = pow(normalizedCenterToRayOrigin, 2) - pow(radius, 2);
    double delta = pow(b, 2) - (4*a*c);

    if (delta < 0) // No intersection, let t be infinity() and return it
        return t1t2;
    else if (delta == 0) // Only 1 intersection
        t1t2.x = (-b) / (2*a);
    else if (delta > 0) { // Two intersections, keep the positive minimum between s1 and s2
        double s1 = (-b - sqrt(delta)) / (2*a);
        double s2 = (-b + sqrt(delta)) / (2*a);
        if (s1 < 0) {
            t1t2.x = s2;
            t1t2.y = s1;
        }
        else {
            t1t2.x = std::min(s1, s2);
            t1t2.y = std::max(s1, s2);
        }
    }

    return t1t2;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Square::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;
    result.t = raySquareIntersection(p0, V);
    result.t2 = std::numeric_limits<double>::infinity();
    if (result.t != std::numeric_limits<double>::infinity()) {
        result.P = p0 + result.t * V;
        result.P2 = 0;
        result.N = normal;
    }
    return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(const vec4& p0, const vec4& V) {
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
            double value1 = dot(normal, cross(v1-v0, P-v0));
            double value2 = dot(normal, cross(v2-v1, P-v1));
            double value3 = dot(normal, cross(v3-v2, P-v2));
            double value4 = dot(normal, cross(v0-v3, P-v3));

            // If all the points are one the same side (all values are positive or negative), then the point is within the rectangle
            if (value1 >= 0 && value2 >= 0 && value3 >= 0 && value4 >= 0)
                return t;
            if (value1 <= 0 && value2 <= 0 && value3 <= 0 && value4 <= 0)
                return t;
        }
    }
    return std::numeric_limits<double>::infinity();
}