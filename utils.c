#include "utils.h"
#include <math.h>

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

Vec4 vec4_sub(Vec4 a, Vec4 b) {
    return (Vec4){a.x - b.x, a.y - b.y, a.z - b.z, 1};
}

Vec4 vec4_add(Vec4 a, Vec4 b) {
    Vec4 result = {a.x + b.x, a.y + b.y, a.z + b.z, 1};
    return result;
}

Vec4 vec4_div(Vec4 a, int v) {
	return (Vec4){a.x / v, a.y / v, a.z / v, 1};
}

Vec4 vec4_mult(Vec4 a, int v) {
	return (Vec4){a.x * v, a.y * v, a.z * v, 1};
}

float vec4_length(Vec4 v)
{
    return sqrtf(v.x * v.x +
                 v.y * v.y +
                 v.z * v.z +
                 v.w * v.w);
}

Vec4 vec4_normalize(Vec4 v)
{
    float len = vec4_length(v);

    if (len > 0.00001f) {
        float inv = 1.0f / len;

        v.x *= inv;
        v.y *= inv;
        v.z *= inv;
        v.w *= inv;
    }
    else {
        v.x = 0.0f;
        v.y = 0.0f;
        v.z = 0.0f;
        v.w = 0.0f;
    }

    return v;
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

float vec3_distance(Vec3 a, Vec3 b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;

    return sqrtf(dx * dx + dy * dy + dz * dz);
}

float vec4_distance(Vec4 a, Vec4 b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;

    return sqrtf(dx * dx + dy * dy + dz * dz);
}

// doesnt use sqrt so faster
float vec4_distance_squared(Vec4 a, Vec4 b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;

    return dx * dx + dy * dy + dz * dz;
}

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
