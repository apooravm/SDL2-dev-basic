#include <SDL2/SDL_error.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_rect.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_timer.h>
#include <SDL2/SDL_video.h>
#include <stdio.h>
#include <SDL2/SDL.h>

#define WIDTH 640
#define HEIGHT 480
#define COLOUR_WHITE 0xffffffff

SDL_Surface* SURFACE;

struct Circle {
	double x;
	double y;
	double r;
};

// homogeneous 4D vector
typedef struct {
	double x, y, z, w;
} Vec4 ;

typedef struct {
	Vec4 vecs[3];
} Triangle;

typedef struct {
	Triangle *tris;
	size_t numTris;
} Mesh;

Mesh* SquareMesh;
Mesh* CubeMesh;

void FillCircle(struct Circle circle, Uint32 colour) {
	double radius_squared = circle.r * circle.r;
	for (double y = circle.y - circle.r; y <= circle.y + circle.r; y++) {
		for (double x = circle.x - circle.r; x <= circle.x + circle.r; x++) {
			double x_x_dist = x - circle.x;
			double y_y_dist = y - circle.y;
			double distance_squared = (x_x_dist * x_x_dist) + (y_y_dist * y_y_dist);

			if (distance_squared < radius_squared) {
				SDL_Rect pixel = (SDL_Rect){x, y, 1, 1};
				SDL_FillRect(SURFACE, &pixel, colour);
			}
		}
	}
}

// convert normalized cartesian point to screen point
void draw_dot(double x, double y) {
	int side = 2;
	// double screenX = (x + 1) * 0.5 * WIDTH;
	// double screenY = (1 - y) * 0.5 * HEIGHT;
	SDL_Rect pixel = (SDL_Rect){x, y, side, side};
	SDL_FillRect(SURFACE, &pixel, COLOUR_WHITE);
}

static inline int norm_to_screen_x(double x) {
    return (int)((x * 0.5 + 0.5) * (WIDTH - 1));
}

static inline int norm_to_screen_y(double y) {
    return (int)(((-y) * 0.5 + 0.5) * (HEIGHT - 1));
}

void DrawLine(double x1, double y1, double x2, double y2) {
    // Convert normalized coords to screen pixels
    int x1i = norm_to_screen_x(x1);
    int y1i = norm_to_screen_y(y1);
    int x2i = norm_to_screen_x(x2);
    int y2i = norm_to_screen_y(y2);

    // Bresenham (unchanged from here)
    int dx = abs(x2i - x1i);
    int dy = abs(y2i - y1i);

    int sx = (x1i < x2i) ? 1 : -1;
    int sy = (y1i < y2i) ? 1 : -1;

    int err = dx - dy;

    while (1) {
        if (x1i >= 0 && x1i < WIDTH && y1i >= 0 && y1i < HEIGHT)
            draw_dot(x1i, y1i);

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
    }
}

int main() {
	SDL_Init(SDL_INIT_VIDEO);
	// window = NULL if SDL_CreateWindow throws error
	SDL_Window* window = SDL_CreateWindow("lmao", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	if (!window) {
		printf("Uh oh\n%s\n", SDL_GetError());
	}

	SDL_Surface* surface = SDL_GetWindowSurface(window);
	SURFACE = surface;

	struct Circle circle = {200, 0, 80};

	int running = 1;
	int move_rate = 1;
	SDL_Event e;

	Vec4 v1 = {0.5, 0, 1, 1};

	while (running) {
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				running = 0;
			}
		}

		double next_py = circle.y + move_rate;
		if (next_py < 480) {
			circle.y = next_py;

		} 

		// draw_dot(200, 200);
		DrawLine(-0.5, 0.5, 0.5, -0.5);
		DrawLine(-1, 0, 1, 0);
		DrawLine(0, 1, 0, -1);

		// FillCircle(surface, circle, COLOUR_WHITE);
		SDL_UpdateWindowSurface(window);
	}

	// SDL_Delay(5000);
	return 0;
}
