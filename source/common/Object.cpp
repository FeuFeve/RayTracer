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
    vec3 rayOrigin = vec3(p0.x, p0.y, p0.z);
    vec3 rayDirection = vec3(V.x, V.y, V.z);
    vec3 centerToRayOrigin = rayOrigin - center;
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

    return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(const vec4& p0, const vec4& V) {
    double t = std::numeric_limits<double>::infinity();
    //TODO: Ray-square intersection;
    return t;
}
