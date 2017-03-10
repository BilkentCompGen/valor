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

#define SONIC_STRAND_FWD 0
#define SONIC_STRAND_REV 1

#define SONIC_GC_WINDOW 100
#define SONIC_GC_SLIDE 100

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
  struct _sonic_interval **gaps;
  struct _sonic_interval **dups;
  struct _sonic_interval **reps;
} sonic;


int  make_sonic(char *, char *, char *, char *, char *);
sonic * load_sonic(char *);
sonic *alloc_sonic(int);
sonic_interval *alloc_sonic_interval(int, int);

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
