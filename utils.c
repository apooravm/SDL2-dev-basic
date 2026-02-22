#include "utils.h"

Vec3 vec3_create(double x, double y, double z) {
    Vec3 v = {x, y, z};
    return v;
}

Vec3 vec3_add(Vec3 a, Vec3 b) {
    Vec3 result = {a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}


// Vec3 vec3_add(Vec3 a, Vec3 b) {
//     Vec3 result = {a.x + b.x, a.y + b.y, a.z + b.z};
//     return result;
// }

// Vec3 vec3_sub(Vec3 a, Vec3 b) {
//     Vec3 result = {a.x - b.x, a.y - b.y, a.z - b.z};
//     return result;
// }

// scale vec3 by s
Vec3 vec3_scale(Vec3 v, float s) { return (Vec3){v.x * s, v.y * s, v.z * s}; }

// Vec3 vec3_scale(Vec3 v, double s) {
//     Vec3 result = {v.x * s, v.y * s, v.z * s};
//     return result;
// }

double vec3_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/* Vec4 functions */

Vec4 vec4_create(double x, double y, double z, double w) {
    Vec4 v = {x, y, z, w};
    return v;
}

Vec4 vec4_from_vec3(Vec3 v, double w) {
    Vec4 result = {v.x, v.y, v.z, w};
    return result;
}

Vec3 vec4_to_vec3(Vec4 v) {
    Vec3 result = {v.x, v.y, v.z};
    return result;
}
