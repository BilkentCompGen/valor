#ifndef __COMMON
#define __COMMON

#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
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

#define MAX_BAMS 256

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
// Track memory usage
extern long long memUsage;
extern FILE *logFile; //Defined in valor.c
extern int CUR_CHR;
typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* reps; /* path to repeatmasker file - *rm.out */
	char* dups; /* path to segmental duplications file - bed */
	char* bam_files; /* paths to comma separated input BAM files as a single string before being tokenized */
	char* bam_list_path; /* path to a file that lists BAM file paths in advance */
	char** bam_file_list; /* the actual list that holds all bam file paths after tokenization */
	char* gaps; /* path to assembly gaps file - bed */
	char* mei;  /* regular expression-like MEI list */
	char* outprefix; /* prefix for the output files */
	int  force_read_length; /* force read length to a certain value, discard those that are shorter. Hidden feature due to GIAB */
	char run_vh; /* boolean stand-in to run VariationHunter */
	char run_rd; /* boolean stand-in to run Read Depth */
	char run_ns; /* boolean stand-in to run NovelSeq */
	char run_sr; /* boolean stand-in to run SPLITREAD */
	char skip_fastq; /* boolean stand-in to skip FASTQ dump */
	char skip_sort; /* boolean stand-in to skip FASTQ sort */
	char skip_remap; /* boolean stand-in to skip FASTQ remap */
	char skip_vhcluster; /* boolean stand-in to skip VH clustering */
	int  threads; /* number of threads to use for parallel mrFAST, and maybe future parallelization of VALOR */
	int num_bams; /* number of input BAM files */
	int quick; /* boolean stand-in to work in bam-only mode (no divet) */
	int ten_x; /*boolean for whether we're using 10x data*/
	char *sonic;
} parameters;

typedef struct _ref_genome
{
	char* ref_name; /* name of the chromosome */
	int chrom_count; /* number of chromosomes */
	int* chrom_lengths; /* lengths of the chromosomes */
	char** chrom_names; /* names of the chromosomes */
	long gen_length; /* total number of base-pairs in the reference genome - total genome length up to chromosome 22 */
	float** gc; /* gc profile for each chromosome */
	int* window_count; /* number of windows for each chromosome */
}ref_genome;





/* Parameter related VALOR functions */
void init_params( parameters**);
void print_params( parameters*);

/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
gzFile safe_fopen_gz( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* General BAM processing functions */
int is_proper( int flag);
int is_discordant( bam1_core_t bam_alignment_core, int min, int max);

int is_alt_concordant( int p1, int p2, int flag, char s1, char s2, int min, int max);
int is_concordant( bam1_core_t bam_alignment_core, int min, int max);
char base_as_char( int base_as_int);
char complement_char( char base);
void qual_to_ascii( char* qual);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_size_int( const void* p, const void* q);

//int count_bed_lines(FILE *);
int findChroIndex(ref_genome* ref, char* chroName);
int isFF_RR(bam1_core_t bam_alignment_core);

// Memory allocation/tracking functions
void* getMem( size_t size);
void resizeMem( void **ptr, size_t old_size, size_t new_size);
void freeMem( void* ptr, size_t size);
double getMemUsage();

#define VALOR_LOG(...) fprintf(logFile,__VA_ARGS__)

int chr_atoi(char *chromosome);
#endif
