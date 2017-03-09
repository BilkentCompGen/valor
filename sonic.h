#ifndef __SONIC
#define __SONIC

#define SONIC_MAGIC 42

#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "sonic_interval.h"

/* Return Codes */
#ifndef RETURN_SUCCESS
#define RETURN_SUCCESS 1
#define RETURN_ERROR 0
#define EXIT_FILE_OPEN_ERROR 5
#define EXIT_SONIC 7
/* Maximum filename length */
#define MAX_LENGTH 1024
#endif

#define STRAND_FWD 0
#define STRAND_REV 1


typedef struct _sonic
{
  int number_of_chromosomes;
  int *number_of_gaps_in_chromosome;
  int *number_of_dups_in_chromosome;
  int *number_of_repeats_in_chromosome;
  int *chromosome_lengths;
  long genome_length;
  char **chromosome_gc_profile;
  char **chromosome_names;
  struct _sonic_interval_array **gaps;
  struct _sonic_interval_array **dups;
  struct _sonic_interval_array **reps;
} sonic;


int  make_sonic(char *, char *, char *, char *, char *);
int load_sonic(char *);

FILE* sonic_fopen( char*, char*);
gzFile sonic_fopen_gz( char*, char*);
int count_bed_lines(FILE *);
int count_bed_chromosome_lines(FILE *, char *);
void sonic_set_str( char**, char*);
void* sonic_get_mem( size_t );
void sonic_write_bed_entries(gzFile, sonic_bed_line *, int, int, char **);
sonic_bed_line *sonic_read_bed_file(FILE *, int, int);
void sonic_write_repeat_item(gzFile, sonic_repeat *);
#endif
