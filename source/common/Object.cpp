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
    return result;
}

/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin center and radius radius    : Find t ------- */
double Sphere::raySphereIntersection(vec4 p0, vec4 V) {
    double t = std::numeric_limits<double>::infinity();
    //TODO: Ray-sphere intersection;

    // t²d.d + 2td.(o-center) + ||o-center||²-r² = 0
    // center = this->center;
    // r = this->radius;
    // o = p0
    // d = V
    // a = d.d
    // b = 2d.(o-center)
    // c = ||o-center||²-r²

    vec3 oMinusC = vec3(p0.x - this->center.x, p0.y - this->center.y, p0.z - this->center.z);

    double a = (V.x * V.x) + (V.y * V.y) + (V.z * V.z) + (V.w * V.w);
    double b = 2 * (V.x * oMinusC.x + V.y * oMinusC.y + V.z * oMinusC.z);
    double c = pow(sqrt(pow(oMinusC.x, 2) + pow(oMinusC.y, 2) + pow(oMinusC.z, 2)), 2) - pow(this->radius, 2);

    double delta = pow(b, 2) - 4 * a * c;
    // if (delta < 0) // No intersection, let t be infinity() and return it
    if (delta == 0) // Only 1 intersection
        t = -b / 2*a;
    else if (delta > 0) { // Two intersections, keep the positive minimum between s1 and s2
        double s1 = (-b - sqrt(delta)) / 2*a;
        double s2 = (-b + sqrt(delta)) / 2*a;
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
double Square::raySquareIntersection(vec4 p0, vec4 V) {
    double t = std::numeric_limits<double>::infinity();
    //TODO: Ray-square intersection;
    return t;
}
