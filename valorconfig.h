#ifndef __VALORCONFIG
#define __VALORCONFIG

#ifndef VALOR_DEFAULT_LOG_FILE
#define VALOR_DEFAULT_LOG_FILE "valor.log"
#endif

#define MAX_FRAG_SIZE  1000 // max segment size (distance between paired end reads)
#define VALOR_FILTER_GAP 1
#define VALOR_FILTER_SAT 1
#define BARCODE_LEN 16
/**************CLONE INFORMATION**************************/
//#define CLONE_MEAN  30000
//#define CLONE_STD_DEV 40000
extern double CLONE_MEAN;
extern double CLONE_STD_DEV;
#define MOLECULE_EXT 80000
#define CLONE_MAX  (CLONE_MEAN + 3 * CLONE_STD_DEV)
#define CLONE_MIN 3 * MAX_FRAG_SIZE // (CLONE_MEAN - 3 * CLONE_STD_DEV)
//#define CLONE_MIN (CLONE_MEAN - 3 * CLONE_STD_DEV)
#define CLONE_MAX_DIST 3000000000
#define CLONE_MIN_DIST 0 
#define MOLECULE_BIN_SIZE 1000
/*************INVERSION INFORMATION****************************/
#define INV_MIN_SIZE  80000 // 80K
#define INV_MAX_SIZE  10000000 // 10M
#define INV_GAP  CLONE_MEAN
#define INV_OVERLAP (-CLONE_MEAN/4) // 1 molecule size
#define INVERSION_MIN_REQUIRED_SUPPORT 8
#define INVERSION_MIN_CLUSTER_SIZE 16
/*************DUPLICATION INFORMION****************************/
#define DUP_OVERLAP (-CLONE_MEAN/4)
#define DUP_GAP CLONE_MEAN
#define DUP_MIN_SIZE CLONE_MEAN //1000
#define DUP_MAX_SIZE 7000000
#define DUP_MAX_DIST 10000000
#define DUP_MIN_DIST 100000
#define VALOR_MOBILE_ELEMENTS "Alu:L1:SVA:HERV"
#define DUPLICATION_MIN_CLUSTER_SIZE 16
#define DUPLICATION_MIN_REQUIRED_SUPPORT 8
/*************DELETION INFORMATION****************************/
#define DELETION_MIN_REQUIRED_SUPPORT 3
#define DELETION_MIN_CLUSTER_SIZE 4
/*************GRAPH PROPERTIES****************************/
#define QCLIQUE_LAMBDA 0.5
#define QCLIQUE_GAMMA 0.6
#define MAX_INVERSIONS_IN_GRAPH 120500
/************** MOLECULE RECOVERY******************/
#define WINDOW_SIZE  (MAX_FRAG_SIZE) // min window size
#define MIN_COVERAGE  0 // min coverage of window size
#define MAX_COVERAGE 50
#define EXTENSION (10 * MAX_FRAG_SIZE) // extension wing
#define MIN_REQUIRED_READS_IN_MOLECULE 8
/*--------------READS---------------------*/
#define MIN_QUAL 0
#define MAX_ALTERNATIVE_CHECK_QUAL 8
#define READ_SAMPLE_SIZE 1000000
#define FILTER1XK 1
#define MAX_SUPPORT 100
#define ALTERNATIVE_MAPPING_BIT 256
#define ALTERNATIVE_MAPPING_FLAG "SA"
#define CHECK_ALTERNATIVE_MAPPINGS 0
void printvalorconfig();

#endif
