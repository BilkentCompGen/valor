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
extern int CUR_CHR;

#define SV_MAX_ID 16
typedef enum SV_TYPE{
        SV_DELETION = 1,
        SV_INVERSION = 2,
        SV_DUPLICATION = 4,
        SV_INVERTED_DUPLICATION = 8,
        SV_TRANSLOCATION = 16,
        SV_INVERTED_TRANSLOCATION = 32
}sv_type;

sv_type atosv(char *str);

typedef struct _params
{
	char* bam_file; /* paths to comma separated input BAM files as a single string before being tokenized */
	char* outprefix; /* prefix for the output files */
	char *sonic_file;
	char *logfile;
	sv_type svs_to_find;
	unsigned int threads; 
	_Bool low_mem;
	int chromosome_count;	 
} parameters;


/* Parameter related VALOR functions */
parameters *get_params(void);
parameters *init_params(void);
void free_params(void *);
void print_params( parameters*);

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

int what_is_min_cluster_size(sv_type type);
int chr_atoi(char *chromosome);
#endif
