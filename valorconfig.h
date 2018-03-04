#ifndef __VALORCONFIG
#define __VALORCONFIG


#define READ_LENGTH  88 // length of each pair end
#define MAX_FRAG_SIZE  1000 // max segment size (distance between paired end reads)
/**************CLONE INFORMATION**************************/
//#define CLONE_MEAN  30000
//#define CLONE_STD_DEV 40000
extern double CLONE_MEAN;
extern double CLONE_STD_DEV;
#define MOLECULE_EXT 90000
#define CLONE_MAX  (CLONE_MEAN + 3 * CLONE_STD_DEV)
#define CLONE_MIN 3 * MAX_FRAG_SIZE // (CLONE_MEAN - 3 * CLONE_STD_DEV)
/*************INVERSION INFORMATION****************************/
#define INV_MIN_SIZE  100000 // 500K
#define INV_MAX_SIZE  10000000 // 10M
#define INV_GAP  CLONE_MEAN
#define INV_OVERLAP (-CLONE_MEAN) // 1 molecule size
/*************DUPLICATION INFORMION****************************/
#define DUP_OVERLAP (-CLONE_MEAN)
#define DUP_GAP CLONE_MEAN
#define DUP_MIN_SIZE 50000
#define DUP_MAX_SIZE 500000
#define DUP_MAX_DIST 10000000
#define DUP_MIN_DIST 50000
/*************GRAPH PROPERTIES****************************/
#define QCLIQUE_LAMBDA 0.35
#define QCLIQUE_GAMMA 0.4
#define QCLIQUE_TABU 20
#define MAX_INVERSIONS_IN_GRAPH 12500
/***************INFER CLONES******************/
#define WINDOW_SIZE  (MAX_FRAG_SIZE) // min window size
#define MIN_COVERAGE  0 // min coverage of window size
#define MAX_COVERAGE 50
#define EXTENSION (10 * MAX_FRAG_SIZE) // extension wing
/*--------------READS---------------------*/
#define MIN_QUAL 0
#define SLIDING_WINDOW  READ_LENGTH
#define READ_SAMPLE_SIZE 1000000
#define FILTER1XK 1
void printvalorconfig();

#endif
