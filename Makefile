dev: main.c
# -lSDL2 linking SDL2; You need to link the SDL2 library explicitly when compiling your program. This can be done by adding -lSDL2 to your gcc command.
	@gcc main.c -o SDL_App_GG.exe -lSDL2 && ./SDL_App_GG.exe

