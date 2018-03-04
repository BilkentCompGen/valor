
#include "progress.h"


void _up_prog(int current, int limit){
        float __ratio;
	if(limit==0){
		__ratio = 1;
	}
	else{
		__ratio = (float)current/limit;
	}
	SETCOLOR(ANSI_YELLOW);
        int __cursor = __ratio * TERMINAL_WIDTH;
        int _i;

        printf("\r");
        for(_i=0;_i<__cursor-1;_i++){
                printf(" ");
        }

        SETCOLOR(ANSI_RED);

        printf(" ");
        SETCOLOR(ANSI_RESET);
        for(;_i<TERMINAL_WIDTH;_i++){
                printf(" ");
        }
        printf(ANSI_UP "\rProgress:%70f%%" ANSI_DOWN,__ratio*100);
        fflush(stdout);

}
