#ifndef __SONIC_REFERENCE
#define __SONIC_REFERENCE

#include "sonic.h"
#include <ctype.h>

int get_number_of_chromosomes(FILE *);
int get_chromosome_info(FILE *, int, int **, char ***);
void sonic_write_gc_profile(gzFile, FILE *, int, char **);
int sonic_find_chromosome_index(char **, const char *, int);
int sonic_refind_chromosome_index(sonic *, const char *);
void  sonic_read_gc_profile(gzFile, sonic *);
#endif
