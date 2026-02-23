#include <stdlib.h>
#include "grid.h"
#include "utils.h"

/* ---------- Internal Helpers ---------- */

static void cell_add_particle(Cell3D* cell, Particle3D* p)
{
    if (cell->count >= cell->capacity) {
        cell->capacity = (cell->capacity == 0) ? 8 : cell->capacity * 2;
        cell->contents = realloc(cell->contents, cell->capacity * sizeof(Particle3D*));
    }

    cell->contents[cell->count++] = p;
}

/* ---------- Grid Creation ---------- */

Grid3D* grid_create(float screenW, float screenH, float screenD, float cellSize)
{
    Grid3D* grid = malloc(sizeof(Grid3D));

    grid->cellSize = cellSize;

    grid->cols   = (int)(screenW / cellSize);
    grid->rows   = (int)(screenH / cellSize);
    grid->layers = (int)(screenD / cellSize);

    grid->cells = malloc(grid->cols * sizeof(Cell3D**));

    for (int i = 0; i < grid->cols; i++) {
        grid->cells[i] = malloc(grid->rows * sizeof(Cell3D*));

        for (int j = 0; j < grid->rows; j++) {
            grid->cells[i][j] = malloc(grid->layers * sizeof(Cell3D));

            for (int k = 0; k < grid->layers; k++) {
                Cell3D* cell = &grid->cells[i][j][k];

                cell->pos.x = i * cellSize;
                cell->pos.y = j * cellSize;
                cell->pos.z = k * cellSize;

                cell->w = cellSize;
                cell->h = cellSize;
                cell->d = cellSize;

                cell->contents = NULL;
                cell->count = 0;
                cell->capacity = 0;
            }
        }
    }

    return grid;
}

void grid_free(Grid3D* grid)
{
    for (int i = 0; i < grid->cols; i++) {
        for (int j = 0; j < grid->rows; j++) {
            for (int k = 0; k < grid->layers; k++) {
                free(grid->cells[i][j][k].contents);
            }

            free(grid->cells[i][j]);
        }

        free(grid->cells[i]);
    }

    free(grid->cells);
    free(grid);
}

/* ---------- Core Operations ---------- */

void grid_reset(Grid3D* grid)
{
    for (int i = 0; i < grid->cols; i++) {
        for (int j = 0; j < grid->rows; j++) {
            for (int k = 0; k < grid->layers; k++) {
                grid->cells[i][j][k].count = 0;
            }
        }
    }
}

void grid_add_particles(Grid3D* grid, Particle3D* particles, int particleCount) {
    for (int n = 0; n < particleCount; n++) {
        Particle3D* p = &particles[n];

        int i = (int)(p->position.x / grid->cellSize);
        int j = (int)(p->position.y / grid->cellSize);
        int k = (int)(p->position.z / grid->cellSize);

        if (i >= 0 && i < grid->cols &&
            j >= 0 && j < grid->rows &&
            k >= 0 && k < grid->layers)
        {
            cell_add_particle(&grid->cells[i][j][k], p);
        }
    }
}

/* ---------- Neighbour Queries ---------- */

Particle3D** grid_get_cell_neighbours(Grid3D* grid, Particle3D* p, int* outCount)
{
    int i = (int)(p->position.x / grid->cellSize);
    int j = (int)(p->position.y / grid->cellSize);
    int k = (int)(p->position.z / grid->cellSize);

    if (i < 0 || i >= grid->cols ||
        j < 0 || j >= grid->rows ||
        k < 0 || k >= grid->layers)
    {
        *outCount = 0;
        return NULL;
    }

    Cell3D* cell = &grid->cells[i][j][k];
    *outCount = cell->count;
    return cell->contents;
}

Particle3D** grid_get_neighbours(Grid3D* grid, Particle3D* p, int* outCount)
{
    int i = (int)(p->position.x / grid->cellSize);
    int j = (int)(p->position.y / grid->cellSize);
    int k = (int)(p->position.z / grid->cellSize);

    int maxPossible = 256;
    Particle3D** neighbours = malloc(maxPossible * sizeof(Particle3D*));
    int count = 0;

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {

                int ni = i + di;
                int nj = j + dj;
                int nk = k + dk;

                if (ni >= 0 && ni < grid->cols &&
                    nj >= 0 && nj < grid->rows &&
                    nk >= 0 && nk < grid->layers)
                {
                    Cell3D* cell = &grid->cells[ni][nj][nk];

                    for (int c = 0; c < cell->count; c++) {
                        Particle3D* p2 = cell->contents[c];

                        if (p2 == p)
                            continue;

						float h2 = grid->cellSize * grid->cellSize;
                        if (vec4_distance_squared(p->position, p2->position) < h2) {
                            neighbours[count++] = p2;
                        }
                    }
                }
            }
        }
    }

    *outCount = count;
    return neighbours;
}
