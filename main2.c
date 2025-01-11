#include <pthread.h>
#include <asm-generic/ioctls.h>
#include <bits/pthreadtypes.h>
#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <time.h>

#define ESC "\033"
#define DO_SOMETHING ESC "[?1049h"  // Switch to alternate screen buffer
#define DO_SOMETHING_OFF ESC "[?1049l"  // Switch back to normal screen buffer
#define MOVE_CURSOR_TO(x, y) ESC "[%d;%dH", x, y
#define CLEAR_SCREEN ESC "[2J"
#define HIDE_CURSOR "\033[?25l"
#define SHOW_CURSOR "\033[?25h"

#define NUM_THREADS 2

int STOP_READING = 1;
pthread_t threads[NUM_THREADS];
pthread_mutex_t mutex;

struct terminal_conf {
	int rows;
	int cols;
	int cursor_pos_x;
	int cursor_pos_y;
	struct termios original_state;
};

struct terminal_conf Term_Conf;

// void* do_nothing() {}

void cleanup_threads() {
    // pthread_mutex_destroy(&mutex);
    // pthread_exit(NULL);
}

void init_thread() {
	pthread_t threads[NUM_THREADS];
	pthread_mutex_init(&mutex, NULL);
}

void die(const char *s) {
	write(STDOUT_FILENO, "\x1b[2J", 4);
	write(STDOUT_FILENO, "\x1b[H", 3);
	perror(s);
	exit(1);
}

void clear_screen() {
	printf(CLEAR_SCREEN);
	fflush(stdout);
}

void waitFor (unsigned int secs) {
    unsigned int retTime = time(0) + secs;   // Get finishing time.
    while (time(0) < retTime);               // Loop until it arrives.
}

// fflush(stdout) is a function in C that forces 
// the output buffer of the stdout stream (standard output) 
// to be written to the terminal (or any other associated output device) immediately.
// move_cursor updates the global x and y cursor positions
void move_cursor(int x, int y) {
	Term_Conf.cursor_pos_x = x;
	Term_Conf.cursor_pos_y = y;

    printf(MOVE_CURSOR_TO(y, x));
	fflush(stdout);
}

// move cursor without updating pos
void move_cursor_NO_REASSGN(int x, int y) {
    printf(MOVE_CURSOR_TO(y, x));
	fflush(stdout);
}

void enable_virtual_window() {
	printf(DO_SOMETHING);
	fflush(stdout);
}

void disable_virtual_window() {
	printf(DO_SOMETHING_OFF);
	fflush(stdout);
}

void disable_raw_mode() {
	cleanup_threads();

	printf(SHOW_CURSOR);
	tcsetattr(STDIN_FILENO, TCSAFLUSH, &Term_Conf.original_state);

    // Switch back to the normal screen buffer
	disable_virtual_window();
}

int get_terminal_size(int *rows, int *cols) {
    struct winsize ws;

	// Catching err
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == -1 || ws.ws_col == 0) {
        return -1;
    }

    *cols = ws.ws_col;
    *rows = ws.ws_row;

	return 0;
}

// Hard way?
int get_terminal_size2(int *rows, int *cols) {
	struct winsize ws;
	if (1 || ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == -1 || ws.ws_col == 0) {
		if (write(STDOUT_FILENO, "\x1b[999C\x1b[999B", 12) != 12) return -1;
		// editorReadKey();
		return -1;
	} else {
		*cols = ws.ws_col;
		*rows = ws.ws_row;
		return 0;
	}
}

// Canonical mode -> raw mode
void enable_raw_mode() {
	enable_virtual_window();
	printf(HIDE_CURSOR);

	tcgetattr(STDIN_FILENO, &Term_Conf.original_state);
	// atexit(disable_raw_mode);

	struct termios raw = Term_Conf.original_state;

	// c_lflag - local flags (misc flags)
	// c_iflag - input flag
	// c_oflag - output flag
	// c_cflag - control flag
	// Turning off
	// ECHO - Typed text wont appear on screen like usual
	// ICANON - Able to read input byte by byte instead of line by line
	raw.c_lflag &= ~(ECHO | ICANON);

	tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
}

void print_live_size(int *p_x, int *p_y) {
	move_cursor(3, 3);
	printf("%d %d", *p_x, *p_y);

	move_cursor(*p_x, *p_y);
}

void print_at_pos(char *text, int pos_x, int pos_y) {
	int temp_x = Term_Conf.cursor_pos_x;
	int temp_y = Term_Conf.cursor_pos_y;

	move_cursor(pos_x, pos_y);
	printf("%s", text);

	// reset to previous pos
	move_cursor(temp_x, temp_y);
}

void init_window() {
	init_thread();
	atexit(disable_raw_mode);

	int rows, cols;
	if (get_terminal_size(&rows, &cols) == -1) {
		die("get_terminal_size");
	}

	Term_Conf.cursor_pos_x = 1;
	Term_Conf.cursor_pos_y = 1;

	Term_Conf.cols = cols;
	Term_Conf.rows = rows;

	move_cursor_NO_REASSGN(1, 3);
	printf("Max cols: %d Max rows: %d ", Term_Conf.cols, Term_Conf.rows);
	move_cursor_NO_REASSGN(Term_Conf.cursor_pos_x, Term_Conf.cursor_pos_y);
}

void draw_rectangle(int x, int y, int l, int b) {
	for (int i_x = 0; i_x < l; i_x++) {
		for (int i_y = 0; i_y <= b; i_y++) {
			pthread_mutex_lock(&mutex);
			move_cursor_NO_REASSGN(x + i_x, y + i_y);
			pthread_mutex_unlock(&mutex);
			printf("#");
		}
	}
}

void* read_user_input(void* thread_id) {
	char c;
	int pos_x = 5;
	int pos_y = 5;

	int rec_x = 4;

	while (read(STDIN_FILENO, &c, 1) == 1) {
		switch (c) {
			case 'q':
				STOP_READING = 0;

			case 'w':
				if (Term_Conf.cursor_pos_y > 1) {
					pos_y--;
				}
				break;
				// printf("W pressed\n");

			case 's':
				if (Term_Conf.cursor_pos_y < Term_Conf.rows) {
					pos_y++;
				}
				break;
				// printf("S pressed\n");
			
			case 'd':
				if (Term_Conf.cursor_pos_x < Term_Conf.cols) {
					pos_x++;
				}
				break;

			case 'a':
				if (Term_Conf.cursor_pos_x > 1) {
					pos_x--;
				}
				break;
		}

		if (!STOP_READING) {
			break;
		}
		// clear_screen();

		pthread_mutex_lock(&mutex);
		move_cursor(pos_x, pos_y);
		pthread_mutex_unlock(&mutex);

		// usleep(100000);

		move_cursor_NO_REASSGN(1, 4);
		printf("X: %d Y: %d", Term_Conf.cursor_pos_x, Term_Conf.cursor_pos_y);
		move_cursor_NO_REASSGN(Term_Conf.cursor_pos_x, Term_Conf.cursor_pos_y);
	}

    pthread_exit(NULL);
}

void* animation(void* thread_id) {
	while (STOP_READING) {
		usleep(100000);
		clear_screen();

		draw_rectangle(Term_Conf.cursor_pos_x, Term_Conf.cursor_pos_y, 4, 9);
	}

    pthread_exit(NULL);
}

int main() {
	enable_raw_mode();
	init_window();

	long read_input_id = 0;
	pthread_create(&threads[read_input_id], NULL, read_user_input, (void*)read_input_id);

	long animation_id = 1;
	pthread_create(&threads[animation_id], NULL, animation, (void*)animation_id);

	pthread_join(threads[animation_id], NULL);
	pthread_join(threads[read_input_id], NULL);

    pthread_mutex_destroy(&mutex); // Clean up the mutex
    pthread_exit(NULL);
}
