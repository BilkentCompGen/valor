#ifndef __COMMON
#define __COMMON
#include "vector.h"
#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
//#include <zlib.h>
#include <omp.h>
#include "valorconfig.h"
#define MAX(a,b) ((a)>(b)?(a):(b)) 
#define MIN(a,b) ((a)<(b)?(a):(b))
/* Exit Codes */
#define EXIT_SUCCESS 0
#define EXIT_COMMON 1
#define EXIT_MAXBAMS 2
#define EXIT_PARAM_ERROR 3
#define EXIT_EXTERNAL_PROG_ERROR 4
#define EXIT_FILE_OPEN_ERROR 5
#define EXIT_READGROUP 6

/* Return Codes */
#define RETURN_SUCCESS 0
#define RETURN_ERROR 1

/* STRAND */
#define READ_STRAND_POS 0 
#define READ_STRAND_NEG 1

/* Maximum filename length */
#define MAX_LENGTH 1024

/* MAPPING INFO */
#define RPUNMAPPED 0
#define RPCONC 1
#define RPDEL 2
#define RPINV 3
#define RPINS 4
#define RPTDUP 5
#define RPMEI 6
#define RPSEND 7
#define RPPP	8
#define RPMM	9
#define RPTDUPPM 10
#define RPTDUPMP 11
#define RPINTER 12


// Track memory usage
extern long long memUsage;
extern FILE *logFile; //Defined in valor.c

#define SV_MAX_ID 32
typedef enum SV_TYPE{
        SV_DELETION = 1,
        SV_TANDEM_DUPLICATION = 2,
        SV_INVERSION = 4,
        SV_DIRECT_DUPLICATION = 8,
        SV_INVERTED_DUPLICATION = 16,
        SV_TRANSLOCATION = 32,
        SV_INVERTED_TRANSLOCATION = 64,
        SV_RECIPROCAL=128,
        SV_INVERTED_RECIPROCAL=256,
}sv_type;


int cmp_int(const void *v1, const void *v2);
sv_type atosv(const char *str);

const char *sv_type_name(sv_type);




/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* General BAM processing functions */
int is_proper( int flag);


int is_alt_concordant( int p1, int p2, int flag, char s1, char s2, int min, int max);
int identify_read_alignment( bam1_core_t bam_alignment_core, int min, int max);
char base_as_char( int base_as_int);
char complement_char( char base);
void qual_to_ascii( char* qual);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_size_int( const void* p, const void* q);




// Memory allocation/tracking functions
void* getMem( size_t size);
void resizeMem( void **ptr, size_t old_size, size_t new_size);
void freeMem( void* ptr, size_t size);
double getMemUsage();

#define VALOR_LOG(...) fprintf(logFile,__VA_ARGS__)

int what_is_min_cluster_size(sv_type type, int ploidy);
int chr_atoi(char *chromosome);
#endif
