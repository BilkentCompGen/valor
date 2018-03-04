#include <stdio.h>
#define ANSI_UP  "\033[1A"
#define ANSI_DOWN "\033[1B"
#define ANSI_RED "\x1b[41m"
#define ANSI_YELLOW "\x1b[43m"
#define ANSI_RESET "\x1b[0m"
#define SETCOLOR(COLOR) do{printf(COLOR);}while(0)
#define TERMINAL_WIDTH 80

void _up_prog(int current, int limit);
#if LIVE_PROGRESS
#define update_progress(c,l) _up_prog(c,l);
#else
#define update_progress(c,l) do{}while(0)
#endif
