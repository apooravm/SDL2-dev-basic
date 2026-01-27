#include <SDL2/SDL.h>
#include <SDL2/SDL_error.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_rect.h>
#include <SDL2/SDL_scancode.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_timer.h>
#include <SDL2/SDL_version.h>
#include <SDL2/SDL_video.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WIDTH 640
#define HEIGHT 480
#define COLOUR_WHITE 0xffffffff
#define COLOUR_RED 0xffff0000
#define COLOUR_GREEN 0x000fff00
#define COLOUR_BLUE 0xff0000ff
#define SQUARE_TRIANGLE_COUNT 2
#define CUBE_TRIANGLE_COUNT 12

SDL_Surface *SURFACE;

struct Circle {
    double x;
    double y;
    double r;
};

typedef struct {
    double x, y, z;
} Vec3;

// homogeneous 4D vector
typedef struct {
    double x, y, z, w;
} Vec4;

typedef struct {
    Vec4 vecs[3];
} Triangle;

typedef struct {
    Triangle *tris;
    size_t numTris;
} Mesh;

// 4x4 matrix
typedef struct {
    double m[4][4];
} Mat4;

typedef struct {
    Vec3 pos;
    double pitch;
    double yaw;

    Vec3 forward;
    Vec3 right;
    Vec3 up;
} Camera;

int D_DIST = 1;
int W_DEF = 1;
double ASPECT_RATIO;
double FOV_ANGLE = 90;
double FOV;
double Z_NORM;
double fNear = 0.1f;
double fFar = 10.0f;
double fFov = 90.0f;
double fFovRad;
Camera camera = {{0.0, 0.0, 1.0}, 0.0, 0.0, 0.0, 0.0, 0.0};
Mat4 projection_matrix = {0};
Mat4 LOOKAT_MTX = {0}; // lookat matrix
float *zbuffer = NULL;

Mesh *SquareMesh;
Mesh *CubeMesh;

void dump_vertex_to_debug_file(float val) {
    FILE *debugFile = fopen("debug.txt", "a");
    if (debugFile == NULL) {
        perror("Failed to open debug file");
        exit(1);
    }

    // Write vertex 0
    fprintf(debugFile, "VAL=%f\n", val);

    fclose(debugFile);
}

void FillCircle(struct Circle circle, Uint32 colour) {
    double radius_squared = circle.r * circle.r;
    for (double y = circle.y - circle.r; y <= circle.y + circle.r; y++) {
        for (double x = circle.x - circle.r; x <= circle.x + circle.r; x++) {
            double x_x_dist = x - circle.x;
            double y_y_dist = y - circle.y;
            double distance_squared =
                (x_x_dist * x_x_dist) + (y_y_dist * y_y_dist);

            if (distance_squared < radius_squared) {
                SDL_Rect pixel = (SDL_Rect){x, y, 1, 1};
                SDL_FillRect(SURFACE, &pixel, colour);
            }
        }
    }
}

Vec4 create_vec4(double x, double y, double z, double w) {
    Vec4 vec = {x, y, z, w};
    return vec;
}
Triangle create_triangle(Vec4 v0, Vec4 v1, Vec4 v2) {
    Triangle tri = {{v0, v1, v2}};
    return tri;
}

// Allocate sq_mesh on its own on heap mem
// Otherwise after the func ends, it frees the obj and pointer is left hanging
Mesh *get_square(double x, double y, double side, double z) {
    Mesh *sq_mesh = (Mesh *)malloc(sizeof(Mesh));
    if (sq_mesh == NULL) {
        printf("Mem allocation failed!\n");
        return NULL;
    }

    // seperate allocation for triangles
    sq_mesh->tris =
        (Triangle *)malloc(SQUARE_TRIANGLE_COUNT * sizeof(Triangle));
    if (sq_mesh->tris == NULL) {
        printf("Memory allocation for triangles failed!\n");
        free(sq_mesh);
        return NULL;
    }

    sq_mesh->numTris = SQUARE_TRIANGLE_COUNT;

    sq_mesh->tris[0] = create_triangle(
        create_vec4(x, y, z, W_DEF), create_vec4(x + side, y + side, z, W_DEF),
        create_vec4(x, y + side, z, W_DEF));

    sq_mesh->tris[1] = create_triangle(
        create_vec4(x, y, z, W_DEF), create_vec4(x + side, y, z, W_DEF),
        create_vec4(x + side, y + side, z, W_DEF));

    return sq_mesh;
}

Mesh *get_cube(double x, double y, double z, double side) {
    // Allocate memory for the cube
    Mesh *cube_mesh = (Mesh *)malloc(sizeof(Mesh));
    if (cube_mesh == NULL) {
        printf("Failed to allocate memory for cube.\n");
        return NULL;
    }

    // Allocate memory for the triangles
    cube_mesh->tris =
        (Triangle *)malloc(CUBE_TRIANGLE_COUNT * sizeof(Triangle));
    if (cube_mesh->tris == NULL) {
        printf("Failed to allocate memory for cube's triangles.\n");
        free(cube_mesh); // Free mesh if triangle allocation fails
        return NULL;
    }

    // Define the number of triangles
    cube_mesh->numTris = CUBE_TRIANGLE_COUNT;

    // Define vertices of the cube
    Vec4 vertices[8] = {
        create_vec4(x, y, z, W_DEF),               // 0: Bottom-front-left
        create_vec4(x + side, y, z, W_DEF),        // 1: Bottom-front-right
        create_vec4(x, y + side, z, W_DEF),        // 2: Top-front-left
        create_vec4(x + side, y + side, z, W_DEF), // 3: Top-front-right
        create_vec4(x, y, z + side, W_DEF),        // 4: Bottom-back-left
        create_vec4(x + side, y, z + side, W_DEF), // 5: Bottom-back-right
        create_vec4(x, y + side, z + side, W_DEF), // 6: Top-back-left
        create_vec4(x + side, y + side, z + side, W_DEF) // 7: Top-back-right
    };

    // Define triangles for each face
    int indices[12][3] = {
        {0, 3, 2}, {0, 1, 3}, // Front face
        {1, 7, 3}, {1, 5, 7}, // Right face
        {5, 6, 7}, {5, 4, 6}, // Back face
        {4, 2, 6}, {4, 0, 2}, // Left face
        {2, 7, 6}, {2, 3, 7}, // Top face
        {4, 1, 0}, {4, 5, 1}  // Bottom face
    };

    for (int i = 0; i < CUBE_TRIANGLE_COUNT; i++) {
        cube_mesh->tris[i] =
            create_triangle(vertices[indices[i][0]], vertices[indices[i][1]],
                            vertices[indices[i][2]]);
    }

    return cube_mesh;
}

// convert normalized cartesian point to screen point
void draw_dot(double x, double y, unsigned int COLOUR) {
    int side = 2;
    // double screenX = (x + 1) * 0.5 * WIDTH;
    // double screenY = (1 - y) * 0.5 * HEIGHT;
    SDL_Rect pixel = (SDL_Rect){x, y, side, side};
    SDL_FillRect(SURFACE, &pixel, COLOUR);
}

static inline int norm_to_screen_x(double x) {
    return (int)((x * 0.5 + 0.5) * (WIDTH - 1));
}

static inline int norm_to_screen_y(double y) {
    return (int)(((-y) * 0.5 + 0.5) * (HEIGHT - 1));
}

int max(int v1, int v2) {
    if (v1 > v2) {
        return v1;
    }

    return v2;
}

void DrawLine(Vec4 *v1, Vec4 *v2) {
    // Convert normalized coords to screen pixels
    int x1i = norm_to_screen_x(v1->x);
    int y1i = norm_to_screen_y(v1->y);
    int x2i = norm_to_screen_x(v2->x);
    int y2i = norm_to_screen_y(v2->y);

    // Bresenham
    int dx = abs(x2i - x1i);
    int dy = abs(y2i - y1i);

    int sx = (x1i < x2i) ? 1 : -1;
    int sy = (y1i < y2i) ? 1 : -1;

    int err = dx - dy;

    double z = v1->z;
    double dz = (v2->z - v1->z) / (double)max(dx, dy);

    while (1) {
        if (x1i >= 0 && x1i < WIDTH && y1i >= 0 && y1i < HEIGHT) {
            int idx = y1i * WIDTH + x1i;
            if (z < zbuffer[idx]) {
                zbuffer[idx] = z;
                draw_dot(x1i, y1i, COLOUR_WHITE);
            }
        }

        if (x1i == x2i && y1i == y2i)
            break;

        int e2 = err * 2;

        if (e2 > -dy) {
            err -= dy;
            x1i += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1i += sy;
        }

        z += dz;
    }
}

void draw_triangle(Triangle *tri) {
    DrawLine(&tri->vecs[0], &tri->vecs[1]);
    DrawLine(&tri->vecs[1], &tri->vecs[2]);
    DrawLine(&tri->vecs[2], &tri->vecs[0]);
}

static inline float edge_function(float ax, float ay, float bx, float by,
                                  float cx, float cy) {
    return (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
}

// great tutorial: https://jtsorlinis.github.io/rendering-tutorial/
void draw_triangle_fill(Triangle *tri, int c_idx) {
    /* Convert NDC -> screen */
    float x0 = (float)norm_to_screen_x(tri->vecs[0].x);
    float y0 = (float)norm_to_screen_y(tri->vecs[0].y);
    float x1 = (float)norm_to_screen_x(tri->vecs[1].x);
    float y1 = (float)norm_to_screen_y(tri->vecs[1].y);
    float x2 = (float)norm_to_screen_x(tri->vecs[2].x);
    float y2 = (float)norm_to_screen_y(tri->vecs[2].y);

    /* Bounding box (floor/ceil, NOT casts) */
    int min_x = (int)floorf(fminf(fminf(x0, x1), x2));
    int max_x = (int)ceilf(fmaxf(fmaxf(x0, x1), x2));
    int min_y = (int)floorf(fminf(fminf(y0, y1), y2));
    int max_y = (int)ceilf(fmaxf(fmaxf(y0, y1), y2));

    /* Clamp to screen */
    if (min_x < 0)
        min_x = 0;
    if (min_y < 0)
        min_y = 0;
    if (max_x >= WIDTH)
        max_x = WIDTH - 1;
    if (max_y >= HEIGHT)
        max_y = HEIGHT - 1;

    /* Triangle area */
    float area = edge_function(x0, y0, x1, y1, x2, y2);
    if (fabsf(area) < 1e-6f)
        return;

    float z0 = tri->vecs[0].z;
    float z1 = tri->vecs[1].z;
    float z2 = tri->vecs[2].z;

    // to correct for the perspective divide done in normalise triangle
    float iz0 = 1.0f / z0;
    float iz1 = 1.0f / z1;
    float iz2 = 1.0f / z2;

    /* Rasterize */
    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {

            float px = (float)x + 0.5f;
            float py = (float)y + 0.5f;

            float w0 = edge_function(x1, y1, x2, y2, px, py);
            float w1 = edge_function(x2, y2, x0, y0, px, py);
            float w2 = edge_function(x0, y0, x1, y1, px, py);

            // barycentric for z values in buffer
            float w0_norm = w0 / area;
            float w1_norm = w1 / area;
            float w2_norm = w2 / area;

            float iz = w0_norm * iz0 + w1_norm * iz1 + w2_norm * iz2;
            float z = 1.0f / iz;
            z = iz;

            float r;
            if (c_idx == 0) {
                r = COLOUR_RED * w0_norm + COLOUR_RED * w1_norm +
                    COLOUR_RED * w2_norm;
            } else {
                r = COLOUR_GREEN * w0_norm + COLOUR_GREEN * w1_norm +
                    COLOUR_GREEN * w2_norm;
            }
            // interpolating depth
            // float z = w0_norm * z0 + w1_norm * z1 + w2_norm * z2;
            // printf("%f\n", tri->vecs[0].x);

            if ((w0 >= 0 && w1 >= 0 && w2 >= 0) ||
                (w0 <= 0 && w1 <= 0 && w2 <= 0)) {
                int idx = y * WIDTH + x;
                if (z < zbuffer[idx]) {
                    zbuffer[idx] = z;
                    draw_dot(x, y, r);
                }
            }
        }
    }
}

// Should NOT update the points in place
// Return the new point to be drawn
Vec4 rotate_xyz(Vec4 *v, double ax, double ay, double az) {
    double x = v->x;
    double y = v->y;
    double z = v->z;

    // Rotate X
    double cosx = cosf(ax);
    double sinx = sinf(ax);
    double y1 = y * cosx - z * sinx;
    double z1 = y * sinx + z * cosx;

    // Rotate Y
    double cosy = cosf(ay);
    double siny = sinf(ay);
    double x2 = x * cosy + z1 * siny;
    double z2 = -x * siny + z1 * cosy;

    // Rotate Z
    double cosz = cosf(az);
    double sinz = sinf(az);
    double x3 = x2 * cosz - y1 * sinz;
    double y3 = x2 * sinz + y1 * cosz;

    Vec4 result;
    result.x = x3;
    result.y = y3;
    result.z = z2;
    result.w = v->w;

    return result;
}

void rotate_triangle(Triangle *tri, double ax, double ay, double az) {
    for (int j = 0; j < 3; j++) {
        tri->vecs[j] = rotate_xyz(&tri->vecs[j], ax, ay, az);
    }
}

Vec4 mat4_mult_vec4_2(const Mat4 *mat, const Vec4 *vec) {
    Vec4 result;

    result.x = vec->x * mat->m[0][0] + vec->y * mat->m[0][1] +
               vec->z * mat->m[0][2] + vec->w * mat->m[0][3];
    result.y = vec->x * mat->m[1][0] + vec->y * mat->m[1][1] +
               vec->z * mat->m[1][2] + vec->w * mat->m[1][3];
    result.z = vec->x * mat->m[2][0] + vec->y * mat->m[2][1] +
               vec->z * mat->m[2][2] + vec->w * mat->m[2][3];
    result.w = vec->x * mat->m[3][0] + vec->y * mat->m[3][1] +
               vec->z * mat->m[3][2] + vec->w * mat->m[3][3];

    // if (result.z <= 0.9f) {
    // 	result.z = 0.9f;
    // }

    return result;
}

void project_triangle(Triangle *tri) {
    tri->vecs[0] = mat4_mult_vec4_2(&projection_matrix, &tri->vecs[0]);
    tri->vecs[1] = mat4_mult_vec4_2(&projection_matrix, &tri->vecs[1]);
    tri->vecs[2] = mat4_mult_vec4_2(&projection_matrix, &tri->vecs[2]);
}

void translate_triangle(Triangle *tri, double val) {
    for (int j = 0; j < 3; j++) {
        tri->vecs[j].z += val;
    }
}

// mat mult like this replaces translating and rotating the triangles separately
void apply_view_matrix(Triangle *tri) {
    tri->vecs[0] = mat4_mult_vec4_2(&LOOKAT_MTX, &tri->vecs[0]);
    tri->vecs[1] = mat4_mult_vec4_2(&LOOKAT_MTX, &tri->vecs[1]);
    tri->vecs[2] = mat4_mult_vec4_2(&LOOKAT_MTX, &tri->vecs[2]);
}

void init_projection_mat() {
    fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * M_PI);
    ASPECT_RATIO = (double)WIDTH / (double)HEIGHT;

    printf("ffovrad: ");
    printf("%f\n", fFovRad);
    printf("Aspect ratio: ");
    printf("%f\n", ASPECT_RATIO);
    printf("Projection [0][0]: ");
    printf("%f\n", fFovRad / ASPECT_RATIO);

    // printf("ffovrad: %f\n", ASPECT_RATIO);
    // ASPECT_RATIO = 1601.0 / 2560.0;
    projection_matrix.m[0][0] = fFovRad / ASPECT_RATIO;
    projection_matrix.m[1][1] = fFovRad;
    projection_matrix.m[2][2] = (fFar + fNear) / (fNear - fFar);
    // projection_matrix.m[2][3] = 1.0f;
    projection_matrix.m[2][3] = (2 * fFar * fNear) / (fNear - fFar);
    // projection_matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
    projection_matrix.m[3][2] = -1.0f;
    projection_matrix.m[3][3] = 0.0f;
}

void normalise_triangle(Triangle *tri) {
    for (int i = 0; i < 3; i++) {
        double w = tri->vecs[i].w;
        if (w != 0.0) {
            tri->vecs[i].x /= w;
            tri->vecs[i].y /= w;
            tri->vecs[i].z /= w;
            tri->vecs[i].w = 1.0;
        }
    }
}

int triangle_in_front(const Triangle *t) {
    for (int i = 0; i < 3; i++) {
        if (t->vecs[i].w > 0)
            return 1;
    }
    return 0;
}

Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

// scale vec3 by s
Vec3 vec3_scale(Vec3 v, float s) { return (Vec3){v.x * s, v.y * s, v.z * s}; }

float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3 cross(Vec3 a, Vec3 b) {
    return (Vec3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x};
}

Vec3 normalize(Vec3 v) {
    float len = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len == 0.0f)
        return v;

    return (Vec3){v.x / len, v.y / len, v.z / len};
}

// Update the View/LookAt matrix every animation loop
// dont understand the math for now :/
void update_camera_basis() {
    Vec3 forward = {
        cosf(camera.pitch) * sinf(camera.yaw), // x
        sinf(camera.pitch),                    // y
        cosf(camera.pitch) * cosf(camera.yaw)  // z
    };
    camera.forward = normalize(forward);

    Vec3 world_up = {0, 1, 0};

    camera.right = normalize(cross(world_up, camera.forward));
    camera.up = cross(camera.forward, camera.right);

    LOOKAT_MTX.m[0][0] = camera.right.x;
    LOOKAT_MTX.m[0][1] = camera.right.y;
    LOOKAT_MTX.m[0][2] = camera.right.z;
    LOOKAT_MTX.m[0][3] = -dot(camera.right, camera.pos);

    LOOKAT_MTX.m[1][0] = camera.up.x;
    LOOKAT_MTX.m[1][1] = camera.up.y;
    LOOKAT_MTX.m[1][2] = camera.up.z;
    LOOKAT_MTX.m[1][3] = -dot(camera.up, camera.pos);

    LOOKAT_MTX.m[2][0] = -camera.forward.x;
    LOOKAT_MTX.m[2][1] = -camera.forward.y;
    LOOKAT_MTX.m[2][2] = -camera.forward.z;
    LOOKAT_MTX.m[2][3] = dot(camera.forward, camera.pos);

    LOOKAT_MTX.m[3][0] = 0;
    LOOKAT_MTX.m[3][1] = 0;
    LOOKAT_MTX.m[3][2] = 0;
    LOOKAT_MTX.m[3][3] = 1;
}

void clear_zbuffer() {
    for (int i = 0; i < WIDTH * HEIGHT; i++)
        zbuffer[i] = 1e9f;
}

static Triangle shipTris[] = {

    /* ===== Body (rectangular hull) ===== */

    // Top
    {{{-0.5, 0.2, -1.0, 1}, {0.5, 0.2, -1.0, 1}, {0.0, 0.2, 1.2, 1}}},
    // Bottom
    {{{0.5, -0.2, -1.0, 1}, {-0.5, -0.2, -1.0, 1}, {0.0, -0.2, 1.2, 1}}},

    // Left side
    {{{-0.5, -0.2, -1.0, 1}, {-0.5, 0.2, -1.0, 1}, {0.0, 0.0, 1.2, 1}}},
    // Right side
    {{{0.5, 0.2, -1.0, 1}, {0.5, -0.2, -1.0, 1}, {0.0, 0.0, 1.2, 1}}},

    /* ===== Nose tip ===== */

    {{{-0.2, 0.1, 1.2, 1}, {0.2, 0.1, 1.2, 1}, {0.0, 0.0, 1.8, 1}}},
    {{{0.2, -0.1, 1.2, 1}, {-0.2, -0.1, 1.2, 1}, {0.0, 0.0, 1.8, 1}}},

    /* ===== Left wing ===== */

    {{{-0.5, 0.0, -0.5, 1}, {-1.4, 0.0, -0.8, 1}, {-0.5, 0.0, 0.4, 1}}},
    {{{-0.5, 0.0, 0.4, 1}, {-1.4, 0.0, -0.8, 1}, {-1.2, 0.0, 0.2, 1}}},

    /* ===== Right wing ===== */

    {{{0.5, 0.0, -0.5, 1}, {0.5, 0.0, 0.4, 1}, {1.4, 0.0, -0.8, 1}}},
    {{{0.5, 0.0, 0.4, 1}, {1.2, 0.0, 0.2, 1}, {1.4, 0.0, -0.8, 1}}},

    /* ===== Engine exhaust ===== */

    {{{-0.3, 0.1, -1.0, 1}, {0.3, 0.1, -1.0, 1}, {0.0, 0.0, -1.6, 1}}},
    {{{0.3, -0.1, -1.0, 1}, {-0.3, -0.1, -1.0, 1}, {0.0, 0.0, -1.6, 1}}},
};

static Mesh shipMesh = {.tris = shipTris,
                        .numTris = sizeof(shipTris) / sizeof(Triangle)};

int main() {
    SDL_Init(SDL_INIT_VIDEO);
    // window = NULL if SDL_CreateWindow throws error
    SDL_Window *window =
        SDL_CreateWindow("lmao", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                         WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (!window) {
        printf("Uh oh\n%s\n", SDL_GetError());
    }

    init_projection_mat();
    zbuffer = malloc(sizeof(float) * WIDTH * HEIGHT);

    SDL_Surface *surface = SDL_GetWindowSurface(window);
    SURFACE = surface;

    // hides cursor, locks it to the window and gives relative motion, xrel and
    // yrel
    SDL_SetRelativeMouseMode(SDL_TRUE);  // FPS mode
    SDL_SetWindowGrab(window, SDL_TRUE); // Optional but recommended
    SDL_ShowCursor(SDL_DISABLE);
    if (SDL_SetRelativeMouseMode(SDL_TRUE) != 0) {
        printf("Relative mode failed: %s\n", SDL_GetError());
    }

    struct Circle circle = {200, 0, 80};

    int running = 1;
    int move_rate = 1;
    SDL_Event e;

    Vec4 v1 = {0.5, 0, 1, 1};

    CubeMesh = get_cube(0, 1.2, 0, 1);
    double angle = 120.0f;

    Mesh *CubeMesh2 = get_cube(0, 0, 0, 1);
    CubeMesh2 = &shipMesh;

    Mesh *Single_tri2 = (Mesh *)malloc(sizeof(Mesh));
    if (Single_tri2 == NULL) {
        printf("Failed to allocate memory for single tri.\n");
        return 1;
    }

    // Allocate memory for the triangles
    Single_tri2->tris = (Triangle *)malloc(1 * sizeof(Triangle));
    if (Single_tri2->tris == NULL) {
        printf("Failed to allocate memory for Single Tri triangles.\n");
        free(Single_tri2); // Free mesh if triangle allocation fails
        return 1;
    }

    Single_tri2->numTris = 1;
    Single_tri2->tris[0] = create_triangle((Vec4){0.0, 0.0, 0.4, W_DEF},
                                           (Vec4){0.0, 1.0, 0.4, W_DEF},
                                           (Vec4){1.0, 0.0, 0.4, W_DEF});

    Mesh *Single_tri = (Mesh *)malloc(sizeof(Mesh));
    if (Single_tri2 == NULL) {
        printf("Failed to allocate memory for single tri.\n");
        return 1;
    }

    // Allocate memory for the triangles
    Single_tri->tris = (Triangle *)malloc(1 * sizeof(Triangle));
    if (Single_tri->tris == NULL) {
        printf("Failed to allocate memory for Single Tri triangles.\n");
        free(Single_tri2); // Free mesh if triangle allocation fails
        return 1;
    }

    Single_tri->numTris = 1;
    Single_tri->tris[0] = create_triangle((Vec4){0.0, 0.0, 0.0, W_DEF},
                                          (Vec4){0.0, 1.0, 0.0, W_DEF},
                                          (Vec4){1.0, 0.0, 0.0, W_DEF});

    int OBJ_SIZE = 1;
    Mesh *objects[OBJ_SIZE];
	objects[0] = CubeMesh;
    // objects[0] = Single_tri2;
    // objects[1] = Single_tri;
    // objects[1] = CubeMesh;
    // objects[1] = CubeMesh2;

    const int FPS = 60;
    const int frameDelay = 1000 / FPS;

    Uint32 frameStart;
    int frameTime;

    Uint64 last = SDL_GetPerformanceCounter();
    float deltaTime = 0.0f;
    float mouseSensitivity = 0.0030f;
    float pitch_yaw_sensitivity = 0.25f;

    while (running) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }

            switch (e.type) {
            case SDL_KEYDOWN:
                if (e.key.keysym.sym == SDLK_ESCAPE) {
                    running = 0;
                }
                break;
            case SDL_QUIT:
                running = 0;
                break;

            case SDL_MOUSEMOTION:
                // break;
                camera.yaw -= e.motion.xrel * mouseSensitivity;
                camera.pitch -= e.motion.yrel * mouseSensitivity;

                if (camera.pitch > 1.55f)
                    camera.pitch = 1.55f;
                if (camera.pitch < -1.55f)
                    camera.pitch = -1.55f;

                break;
            }
        }

        SDL_FillRect(surface, NULL, SDL_MapRGB(surface->format, 0, 0, 0));

        frameStart = SDL_GetTicks();

        const Uint8 *keys = SDL_GetKeyboardState(NULL);

        Uint64 now = SDL_GetPerformanceCounter();
        deltaTime = (float)(now - last) / SDL_GetPerformanceFrequency();
        last = now;

        update_camera_basis();

        float moveSpeed = 5.0f * deltaTime;

        if (keys[SDL_SCANCODE_W]) {
            camera.pos =
                vec3_sub(camera.pos, vec3_scale(camera.forward, moveSpeed));
        }
        if (keys[SDL_SCANCODE_S]) {
            camera.pos =
                vec3_add(camera.pos, vec3_scale(camera.forward, moveSpeed));
        }
        if (keys[SDL_SCANCODE_A]) {
            camera.pos =
                vec3_add(camera.pos, vec3_scale(camera.right, moveSpeed));
        }
        if (keys[SDL_SCANCODE_D]) {
            camera.pos =
                vec3_sub(camera.pos, vec3_scale(camera.right, moveSpeed));
        }
        if (keys[SDL_SCANCODE_L]) {
            camera.yaw += moveSpeed * pitch_yaw_sensitivity;
        }
        if (keys[SDL_SCANCODE_K]) {
            // camera.yaw += moveSpeed * pitch_yaw_sensitivity;
            camera.pitch -= moveSpeed * pitch_yaw_sensitivity;

            if (camera.pitch > 1.55f)
                camera.pitch = 1.55f;
            if (camera.pitch < -1.55f)
                camera.pitch = -1.55f;
        }
        if (keys[SDL_SCANCODE_H]) {
            camera.yaw -= moveSpeed * pitch_yaw_sensitivity;
        }
        if (keys[SDL_SCANCODE_J]) {
            // camera.yaw += moveSpeed * pitch_yaw_sensitivity;
            camera.pitch += moveSpeed * pitch_yaw_sensitivity;

            if (camera.pitch > 1.55f)
                camera.pitch = 1.55f;
            if (camera.pitch < -1.55f)
                camera.pitch = -1.55f;
        }

        double next_py = circle.y + move_rate;
        if (next_py < 480) {
            circle.y = next_py;
        }

        // draw polygons
        for (int m = 0; m < OBJ_SIZE; m++) {
            Mesh *mesh = objects[m];
            for (int i = 0; i < mesh->numTris; i++) {
                Triangle tri_updated = mesh->tris[i];
                // rotate_triangle(&tri_updated, 0, 0.0 * 0.5, 0 * 0.33);
                // translate_triangle(&tri_updated, 1.0);
                apply_view_matrix(&tri_updated);
                project_triangle(&tri_updated);
                if (triangle_in_front(&tri_updated)) {
                    continue;
                }
                normalise_triangle(&tri_updated);
                draw_triangle_fill(&tri_updated, m); // rasterization
            }
        }

        clear_zbuffer();

        angle += 0.02;

        // FillCircle(surface, circle, COLOUR_WHITE);
        SDL_UpdateWindowSurface(window);
        frameTime = SDL_GetTicks() - frameStart;

        if (frameDelay > frameTime) {
            SDL_Delay(frameDelay - frameTime);
        }
    }

    // free(CubeMesh2->tris);
    // free(CubeMesh2);

    free(CubeMesh->tris);
    free(CubeMesh);

    free(zbuffer);
	free(Single_tri->tris);
	free(Single_tri);
	free(Single_tri2->tris);
	free(Single_tri2);

    return 0;
}

/**
      case SDL_KEYDOWN:
        if (e.key.keysym.sym == SDLK_ESCAPE) {
          running = 0;
        }

        if (e.key.keysym.sym == SDLK_w) {
          printf("w key pressed\n");
        }
        break;

      case SDL_KEYUP:
        if (e.key.keysym.sym == SDLK_w) {
          printf("w key released\n");
        }
        break;

      case SDL_MOUSEMOTION:
        printf("Mouse moved: x=%d y=%d | dx=%d dy=%d\n", e.motion.x, e.motion.y,
               e.motion.xrel, e.motion.yrel);
        break;

      case SDL_MOUSEBUTTONDOWN:
        if (e.button.button == SDL_BUTTON_LEFT) {
          printf("Left mouse button pressed at %d, %d\n", e.button.x,
                 e.button.y);
        }
        break;

      case SDL_MOUSEBUTTONUP:
        if (e.button.button == SDL_BUTTON_LEFT) {
          printf("Left mouse button released\n");
        }
        break;
        **/
