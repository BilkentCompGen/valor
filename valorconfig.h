#ifndef __VALORCONFIG
#define __VALORCONFIG
#ifndef VALOR_DEFAULT_LOG_FILE
#define VALOR_DEFAULT_LOG_FILE "valor.log"
#endif
#include <stdio.h>

void printvalorconfig(FILE *file);
extern double SV_OVERLAP_RATIO;
#define CLONE_MAX   (CLONE_MEAN + 3 * CLONE_STD_DEV)
#define CLONE_MIN   (3 * MAX_FRAG_SIZE)
extern int MAX_FRAG_SIZE;
extern int VALOR_FILTER_GAP;
extern int VALOR_FILTER_SAT;
extern int BARCODE_LEN;
/**************CLONE INFORMATION**************************/
extern double CLONE_MEAN;
extern double CLONE_STD_DEV;
extern int MOLECULE_EXT;

extern long CLONE_MAX_DIST;
extern int CLONE_MIN_DIST;
extern int MOLECULE_BIN_SIZE;
/*************INVERSION INFORMATION****************************/
extern int INV_MIN_SIZE;
extern int INV_MAX_SIZE;
#define INV_GAP (SV_OVERLAP_RATIO * CLONE_MEAN)
#define INV_OVERLAP (- SV_OVERLAP_RATIO * CLONE_MEAN)
extern int INVERSION_MIN_REQUIRED_SUPPORT;
extern int INVERSION_MIN_CLUSTER_SIZE;
/*************DUPLICATION INFORMION****************************/

#define DUP_GAP (SV_OVERLAP_RATIO * CLONE_MEAN)
#define DUP_OVERLAP (- SV_OVERLAP_RATIO * CLONE_MEAN)

extern int DUP_MIN_SIZE;
extern int DUP_MAX_SIZE;
extern int DUP_MAX_DIST;
extern int DUP_MIN_DIST;
/*************INTER TRANSLOCATION INFORMION****************************/
#define TRA_GAP (SV_OVERLAP_RATIO * CLONE_MEAN)
#define TRA_OVERLAP (- SV_OVERLAP_RATIO * CLONE_MEAN)
extern int TRA_MIN_SIZE;
extern int TRA_MAX_SIZE;
extern char* VALOR_MOBILE_ELEMENTS;

extern int DUPLICATION_MIN_CLUSTER_SIZE;
extern int TRANSLOCATION_MIN_CLUSTER_SIZE;
extern int DUPLICATION_MIN_REQUIRED_SUPPORT;

extern int TANDEM_DUPLICATION_MIN_CLUSTER_SIZE;
extern int TANDEM_DUPLICATION_MIN_SUPPORT;
/*************DELETION INFORMATION****************************/
extern int DELETION_MIN_REQUIRED_SUPPORT;
extern int DELETION_MIN_CLUSTER_SIZE;

/*************INTER_CHR_EVENTS****************************/
extern int MIN_INTER_CLUSTER_SIZE;
extern int TRA_MIN_INTRA_SPLIT;
/*************GRAPH PROPERTIES****************************/
extern double QCLIQUE_LAMBDA;
extern double QCLIQUE_GAMMA;
#define MAX_INVERSIONS_IN_GRAPH 120500
/************** MOLECULE RECOVERY******************/

extern int MIN_COVERAGE;
extern int MAX_COVERAGE;
extern int EXTENSION;
extern int MIN_REQUIRED_READS_IN_MOLECULE;
/*--------------READS---------------------*/
extern int MIN_QUAL;
extern int MAX_ALTERNATIVE_CHECK_QUAL;
extern int READ_SAMPLE_SIZE;
extern int FILTER1XK;
extern int MAX_SUPPORT;
extern int ALTERNATIVE_MAPPING_BIT;
extern char* ALTERNATIVE_MAPPING_FLAG;
extern int CHECK_ALTERNATIVE_MAPPINGS;

#endif

