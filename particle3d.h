#ifndef PARTICLE3D_H
#define PARTICLE3D_H

#include "utils.h"
#include <SDL2/SDL.h>


/* =========================
   Particle3D Structure
   ========================= */
typedef struct {

    /* Core physics vectors */
    Vec4 position;       // w = 1.0
    Vec4 velocity;       // w = 0.0
    Vec4 acceleration;   // w = 0.0

    /* Physical properties */
    double mass;
    double density;

    /* Rendering / limits */
    double size;
    double maxVel;
    double frictionRate;

    /* Screen bounds */
    int screenW;
    int screenH;
    int screenD;

    /* Containing cube */
    int cubeX;
    int cubeY;
    int cubeZ;
    int cubeWidth;

} Particle3D;


/* =========================
   Function Prototypes
   ========================= */

void Particle3D_Init(Particle3D *p,
                     double x, double y, double z,
                     double size,
                     int screenW, int screenH, int screenD,
                     int cubeX, int cubeY, int cubeZ,
                     int cubeWidth);

void Particle3D_AddForce(Particle3D *p,
                         double fx,
                         double fy,
                         double fz);

void Particle3D_Update(Particle3D *p,
                       double dt);

void Particle3D_Render(Particle3D *p,
                       SDL_Renderer *renderer);

#endif
