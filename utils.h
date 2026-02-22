#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x, y, z;
} Vec3;

// homogeneous 4D vector
typedef struct {
    double x, y, z, w;
} Vec4;

/* Vec3 functions */
Vec3 vec3_create(double x, double y, double z);
Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_sub(Vec3 a, Vec3 b);
Vec3 vec3_scale(Vec3 v, float s);
double vec3_dot(Vec3 a, Vec3 b);

/* Vec4 functions */
Vec4 vec4_create(double x, double y, double z, double w);
Vec4 vec4_from_vec3(Vec3 v, double w);
Vec3 vec4_to_vec3(Vec4 v);

#ifdef __cplusplus
}
#endif

#endif // UTILS_H
