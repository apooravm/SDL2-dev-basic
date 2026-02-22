#ifndef PARTICLE_H
#define PARTICLE_H

#include <SDL2/SDL.h>
typedef struct {
    float x, y;   // position
    float vx, vy; // velocity
    float ax, ay; // acceleration
    float mass;
    float density;
    float width; // square size
    float maxVel;
    float frictionRate;
    int screenW;
    int screenH;
} Particle;

#endif

void Particle_Init(Particle *p, float x, float y, float width, int screenW,
                   int screenH) {
    p->x = x;
    p->y = y;
    p->vx = 0.0f;
    p->vy = 0.0f;
    p->ax = 0.0f;
    p->ay = 0.0f;
    p->mass = 1.0f;
    p->density = 1.0f;
    p->width = width;
    p->maxVel = 40.0f;
    p->frictionRate = 1.0f;
    p->screenW = screenW;
    p->screenH = screenH;
}

void Particle_AddForce(Particle *p, float fx, float fy) {
    p->ax += fx;
    p->ay += fy;
}

void Particle_Update(Particle *p, float dt) {

    // Update velocity
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;

    // Clamp velocity
    if (p->vx > p->maxVel)
        p->vx = p->maxVel;
    if (p->vx < -p->maxVel)
        p->vx = -p->maxVel;
    if (p->vy > p->maxVel)
        p->vy = p->maxVel;
    if (p->vy < -p->maxVel)
        p->vy = -p->maxVel;

    // Apply friction
    p->vx *= p->frictionRate;
    p->vy *= p->frictionRate;

    // Update position
    p->x += p->vx * dt;
    p->y += p->vy * dt;

    float half = p->width / 2.0f;

    // Bounce X
    if (p->x > p->screenW - half) {
        p->x = p->screenW - half;
        p->vx *= -p->frictionRate;
    } else if (p->x < half) {
        p->x = half;
        p->vx *= -p->frictionRate;
    }

    // Bounce Y
    if (p->y > p->screenH - half) {
        p->y = p->screenH - half;
        p->vy *= -p->frictionRate;
    } else if (p->y < half) {
        p->y = half;
        p->vy *= -p->frictionRate;
    }

    // Reset acceleration
    p->ax = 0.0f;
    p->ay = 0.0f;
}

void Particle_Render(Particle *p, SDL_Renderer *renderer) {
    SDL_Rect rect;
    rect.w = (int)p->width;
    rect.h = (int)p->width;
    rect.x = (int)(p->x - p->width / 2.0f);
    rect.y = (int)(p->y - p->width / 2.0f);

    SDL_RenderFillRect(renderer, &rect);
}

