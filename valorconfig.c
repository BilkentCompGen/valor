

#include <stdio.h>
#include "valorconfig.h"

void printvalorconfig(){
	printf("WINDOW_SIZE : %d\n" , WINDOW_SIZE);
	printf("MIN_COVERAGE : %f\n" , MIN_COVERAGE);
	printf("EXTENSION : %d\n" , EXTENSION);

        printf("**************** RUN  DETAILS ********************");
		printf("READ_LENGTH : %d\n" , READ_LENGTH);
        printf("FRAG_SIZE : %d\n" , FRAG_SIZE);

        //printf("CLONE_SIZE : %d\n" , CLONE_SIZE);
        printf("CLONE_MEAN : %d\n" , CLONE_MEAN);
        printf("CLONE_STD_DEV : %d\n" , CLONE_STD_DEV);
        printf("CLONE_MAX : %d\n" , CLONE_MAX);
        printf("CLONE_MIN : %d\n" , CLONE_MIN);

        printf("INV_MIN_SIZE : %d\n" , INV_MIN_SIZE);
        printf("INV_MAX_SIZE : %d\n" , INV_MAX_SIZE);
        printf("INV_GAP : %d\n" , INV_GAP);
        printf("INV_OVERLAP : %d\n" , INV_OVERLAP);

        printf("QCLIQUE_LAMBDA : %f\n" , QCLIQUE_LAMBDA);
	printf("QCLIQUE_GAMMA : %f\n" , QCLIQUE_GAMMA);
	printf("*************************************************"); 
}

