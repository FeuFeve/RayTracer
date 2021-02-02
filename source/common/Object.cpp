//////////////////////////////////////////////////////////////////////////////
//
//  --- Object.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////


#include "common.h"

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Sphere::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;

    //TODO: Ray-sphere setup
    result.t = raySphereIntersection(p0, V);
    if (result.t == std::numeric_limits<double>::infinity())
        return result;

    // r(t) = o + td
    // P = r(t), o = p0, d = V
    result.P = p0 + (result.t * V);
    result.N = p0 - center;

    return result;
}

/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin center and radius radius    : Find t ------- */
double Sphere::raySphereIntersection(const vec4& p0, const vec4& V) {
    double t = std::numeric_limits<double>::infinity();

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
        return t;
    else if (delta == 0) // Only 1 intersection
        t = (-b) / (2*a);
    else if (delta > 0) { // Two intersections, keep the positive minimum between s1 and s2
        double s1 = (-b - sqrt(delta)) / (2*a);
        double s2 = (-b + sqrt(delta)) / (2*a);
        if (s1 < 0)
            t = s2;
        else
            t = std::min(s1, s2);
    }

    return t;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Square::intersect(vec4 p0, vec4 V) {
    IntersectionValues result;

    //TODO: Ray-square setup
    result.t = raySquareIntersection(p0, V);
    if (result.t == std::numeric_limits<double>::infinity())
        return result;

    // r(t) = o + td
    // P = r(t), o = p0, d = V
    result.P = p0 + (result.t * V);
    result.N = normal;

    return result;

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(const vec4& p0, const vec4& V) {
//    double t = std::numeric_limits<double>::infinity();

    vec3 a = vec3(mesh.vertices[5].x, mesh.vertices[5].y, mesh.vertices[5].z); //a.w = 0;// up left
    vec3 b = vec3(mesh.vertices[1].x, mesh.vertices[1].y, mesh.vertices[1].z); //b.w = 0;// up right
    vec3 c = vec3(mesh.vertices[2].x, mesh.vertices[2].y, mesh.vertices[2].z); //c.w = 0;// down right
    vec3 planeNormal = cross(b - a, c - a);
//    printf("planeNormal=(%f, %f, %f)\n", planeNormal.x, planeNormal.y, planeNormal.z);

    //TODO: Ray-square intersection;
    vec3 rayOrigin = vec3(p0.x, p0.y, p0.z); // o
    vec3 rayDirection = vec3(V.x, V.y, V.z); // d
    double distance = point.x*normal.x + point.y*normal.y + point.z*normal.z; // D

    double denom = dot(planeNormal, rayDirection);
    double denom2 = dot(normal, rayDirection);
    printf("denom=%f, denom2=%f\n", denom, denom2);

    double dDotN = dot(rayDirection, normal); // d.n
    if (dDotN == 0) // The ray is parallel and distinct from the plane
    {
        printf("t == +infinity (parallel)\n\n");
        return std::numeric_limits<double>::infinity();
    }

    double t = (distance * dot(rayOrigin, normal)) / dDotN;
//    printf("o=(%f, %f, %f), d=(%f, %f, %f), a=(%f, %f, %f), n=(%f, %f, %f)\n", p0.x, p0.y, p0.z, V.x, V.y, V.z, point.x, point.y, point.z, normal.x, normal.y, normal.z);
//    printf("distance=%f, dot(rayOrigin, normal)=%f, dDotN=%f\n", distance, dot(rayOrigin, normal), dDotN);
    if (t > 0) // Intersection in front of the camera
    {
        printf("t > 0 (intersection in front of the camera)\n\n");
        return t;
    }
    else // Intersection behind the camera, or no intersection
    {
        if (t == 0)
            printf("t == 0 (ray into plane)\n\n");
        else
            printf("t < 0 (intersection behind the camera\n\n");
        return std::numeric_limits<double>::infinity();
    }
}
