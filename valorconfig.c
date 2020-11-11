

#include <stdio.h>
#include "valorconfig.h"


#define GET_FIELD(F) #F
void printvalorconfig(FILE *file){

    fprintf(file,"**************** RUN  DETAILS ********************");

    fprintf(file,"FILTER GAP: %d\n" , VALOR_FILTER_GAP);
    fprintf(file,"FILTER SATELLITE: %d\n" , VALOR_FILTER_SAT);
    fprintf(file,"BARCODE LEN: %d\n" , BARCODE_LEN);
    fprintf(file,"SPLIT MOLECULE MAX DISTANCE: %ld\n" , CLONE_MAX_DIST);
    fprintf(file,"SPLIT MOLECULE MIN DISTANCE: %d\n" , CLONE_MIN_DIST);
    fprintf(file,"MOLECULE DEPT BIN SIZE: %d\n" , MOLECULE_BIN_SIZE);


    fprintf(file,"MIN COVERAGE : %d\n" , MIN_COVERAGE);
    fprintf(file,"MAX COVERAGE : %d\n" , MAX_COVERAGE);

    fprintf(file,"EXTENSION : %d\n" , EXTENSION);
    fprintf(file,"MIN REQUIRED READS FOR A MOLECULE : %d\n" , MIN_REQUIRED_READS_IN_MOLECULE);


    fprintf(file,"MIN QUAL: %d\n" , MIN_QUAL);

    fprintf(file,"Stat sampling size: %d\n" , READ_SAMPLE_SIZE);

    if(FILTER1XK){
        fprintf(file,"Max Svs in the graph: %d\n" , MAX_INVERSIONS_IN_GRAPH);
    }
   fprintf(file,"Max discordant support: %d\n", MAX_SUPPORT); 



    fprintf(file,"MAX_FRAG_SIZE : %d\n" , MAX_FRAG_SIZE);




    fprintf(file,"DUP_MAX_SIZE : %d\n" , DUP_MAX_SIZE);
    fprintf(file,"DUP_MIN_SIZE : "  GET_FIELD(DUP_MIN_SIZE) "\n");
    fprintf(file,"DUP_MAX_DIST : %d\n" , DUP_MAX_DIST);
    fprintf(file,"DUP_MIN_DIST : %d\n" , DUP_MIN_DIST);
    

    fprintf(file,"Duplication minimum cluster size : %d\n" , DUPLICATION_MIN_CLUSTER_SIZE);
    fprintf(file,"Duplication minimum required support : %d\n" , DUPLICATION_MIN_REQUIRED_SUPPORT);
    fprintf(file,"DUP_MIN_DIST : %d\n" , DUP_MIN_DIST);

    fprintf(file,"INV_MIN_SIZE : %d\n" , INV_MIN_SIZE);
    fprintf(file,"INV_MAX_SIZE : %d\n" , INV_MAX_SIZE);


    fprintf(file,"Inversion minimum cluster size : %d\n" , INVERSION_MIN_CLUSTER_SIZE);
    fprintf(file,"Inversion minimum required support : %d\n" , INVERSION_MIN_REQUIRED_SUPPORT);
    
    fprintf(file,"Deletion minimum cluster size : %d\n" , DELETION_MIN_CLUSTER_SIZE);
    fprintf(file,"Deletion minimum required support : %d\n" , DELETION_MIN_REQUIRED_SUPPORT);

    
    fprintf(file,"QCLIQUE_LAMBDA : %f\n" , QCLIQUE_LAMBDA);
    fprintf(file,"QCLIQUE_GAMMA : %f\n" , QCLIQUE_GAMMA);
    fprintf(file,"*************************************************"); 
}

