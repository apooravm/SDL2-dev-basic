#include "particle3d.h"
#include <stdlib.h>

void Particle3D_Init(Particle3D *p, double x, double y, double z, double size,
                     int screenW, int screenH, int screenD, int cubeX,
                     int cubeY, int cubeZ, int cubeWidth) {
	int rand_X = cubeX + rand() % (cubeX + cubeWidth - cubeX + 1);
	int rand_Y = cubeY + rand() % (cubeY + cubeWidth - cubeY + 1);
	int rand_Z = cubeZ + rand() % (cubeZ + cubeWidth - cubeZ + 1);

    p->position = (Vec4){rand_X, rand_Y, rand_Z, 1.0};

    double mult_fact = 1;
    p->velocity = (Vec4){((double)rand() / RAND_MAX - 0.5) * mult_fact,
                         ((double)rand() / RAND_MAX - 0.5) * mult_fact,
                         ((double)rand() / RAND_MAX - 0.5) * mult_fact, 0.0};

    p->acceleration = (Vec4){0.0, 0.0, 0.0, 0.0};

    p->mass = 1.0;
    p->density = 1.0;

    p->size = size;
    p->maxVel = 40.0;
    p->frictionRate = 1.0;

    p->screenW = screenW;
    p->screenH = screenH;
    p->screenD = screenD;

    p->cubeX = cubeX;
    p->cubeY = cubeY;
    p->cubeZ = cubeZ;
    p->cubeWidth = cubeWidth;
}

void Particle3D_AddForce(Particle3D *p, double fx, double fy, double fz) {
    p->acceleration.x += fx / p->mass;
    p->acceleration.y += fy / p->mass;
    p->acceleration.z += fz / p->mass;
}

void Particle3D_Update(Particle3D *p, double dt) {
    // Update velocity
    p->velocity.x += p->acceleration.x * dt;
    p->velocity.y += p->acceleration.y * dt;
    p->velocity.z += p->acceleration.z * dt;

    // Clamp velocity
    if (p->velocity.x > p->maxVel)
        p->velocity.x = p->maxVel;
    if (p->velocity.x < -p->maxVel)
        p->velocity.x = -p->maxVel;

    if (p->velocity.y > p->maxVel)
        p->velocity.y = p->maxVel;
    if (p->velocity.y < -p->maxVel)
        p->velocity.y = -p->maxVel;

    if (p->velocity.z > p->maxVel)
        p->velocity.z = p->maxVel;
    if (p->velocity.z < -p->maxVel)
        p->velocity.z = -p->maxVel;

    // Apply friction
    p->velocity.x *= p->frictionRate;
    p->velocity.y *= p->frictionRate;
    p->velocity.z *= p->frictionRate;

    // Update position
    p->position.x += p->velocity.x * dt;
    p->position.y += p->velocity.y * dt;
    p->position.z += p->velocity.z * dt;

    double half = p->size / 2.0;

    // Bounce X
    if (p->position.x > p->cubeX + p->cubeWidth - half) {
        p->position.x = p->cubeX + p->cubeWidth - half;
        p->velocity.x *= -p->frictionRate;
    } else if (p->position.x < p->cubeX + half) {
        p->position.x = p->cubeX + half;
        p->velocity.x *= -p->frictionRate;
    }

    // Bounce Y
    if (p->position.y > p->cubeY + p->cubeWidth - half) {
        p->position.y = p->cubeY + p->cubeWidth - half;
        p->velocity.y *= -p->frictionRate;
    } else if (p->position.y < p->cubeY + half) {
        p->position.y = p->cubeY + half;
        p->velocity.y *= -p->frictionRate;
    }

    // Bounce Z
    if (p->position.z > p->cubeZ + p->cubeWidth - half) {
        p->position.z = p->cubeZ + p->cubeWidth - half;
        p->velocity.z *= -p->frictionRate;
    } else if (p->position.z < p->cubeZ + half) {
        p->position.z = p->cubeZ + half;
        p->velocity.z *= -p->frictionRate;
    }

    // Reset acceleration
    p->acceleration.x = 0.0;
    p->acceleration.y = 0.0;
    p->acceleration.z = 0.0;
}

void Particle3D_Render(Particle3D *p, SDL_Renderer *renderer) {
    double cameraDistance = 500.0;
    double perspective = cameraDistance / (cameraDistance + p->position.z);

    double projectedX = p->position.x * perspective;
    double projectedY = p->position.y * perspective;
    double projectedSize = p->size * perspective;

    SDL_Rect rect;
    rect.w = (int)projectedSize;
    rect.h = (int)projectedSize;
    rect.x = (int)(projectedX - projectedSize / 2.0);
    rect.y = (int)(projectedY - projectedSize / 2.0);

    SDL_RenderFillRect(renderer, &rect);
}
