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

struct Circle {
	double x;
	double y;
	double r;
};

// homogeneous 4D vector
typedef struct {
	double x, y, z, w;
} Vec4 ;

void FillCircle(SDL_Surface* surface, struct Circle circle, Uint32 colour) {
	double radius_squared = circle.r * circle.r;
	for (double y = circle.y - circle.r; y <= circle.y + circle.r; y++) {
		for (double x = circle.x - circle.r; x <= circle.x + circle.r; x++) {
			double x_x_dist = x - circle.x;
			double y_y_dist = y - circle.y;
			double distance_squared = (x_x_dist * x_x_dist) + (y_y_dist * y_y_dist);

			if (distance_squared < radius_squared) {
				SDL_Rect pixel = (SDL_Rect){x, y, 1, 1};
				SDL_FillRect(surface, &pixel, colour);
			}
		}
	}
}

void draw_dot(const Vec4* v, SDL_Surface* surface) {
	int side = 5;
	double screenX = (v->x + 1) * 0.5 * WIDTH;
	double screenY = (1 - v->y) * 0.5 * HEIGHT;
	SDL_Rect pixel = (SDL_Rect){screenX, screenY, side, side};
	SDL_FillRect(surface, &pixel, COLOUR_WHITE);
}

int main() {
	SDL_Init(SDL_INIT_VIDEO);
	// window = NULL if SDL_CreateWindow throws error
	SDL_Window* window = SDL_CreateWindow("lmao", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	if (!window) {
		printf("Uh oh\n%s\n", SDL_GetError());
	}

	SDL_Surface* surface = SDL_GetWindowSurface(window);

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

		draw_dot(&v1, surface);

		// FillCircle(surface, circle, COLOUR_WHITE);
		SDL_UpdateWindowSurface(window);
	}

	// SDL_Delay(5000);
	return 0;
}
