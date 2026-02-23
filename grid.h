#ifndef GRID_H
#define GRID_H

#include "particle3d.h"

typedef struct {
    Vec3 pos;
    float w, h, d;

    Particle3D** contents;
    int count;
    int capacity;

} Cell3D;

typedef struct {
    Cell3D*** cells;

    int cols;
    int rows;
    int layers;

    float cellSize;   // smoothing radius
} Grid3D;


/* Creation / Destruction */
Grid3D* grid_create(float screenW, float screenH, float screenD, float cellSize);
void grid_free(Grid3D* grid);

/* Core operations */
void grid_reset(Grid3D* grid);
void grid_add_particles(Grid3D* grid, Particle3D* particles, int particleCount);

/* Neighbour queries */
Particle3D** grid_get_cell_neighbours(Grid3D* grid, Particle3D* p, int* outCount);
Particle3D** grid_get_neighbours(Grid3D* grid, Particle3D* p, int* outCount);

#endif
