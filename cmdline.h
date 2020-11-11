#ifndef CMDLINE_H
#define CMDLINE_H
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vector.h"
#include "valorconfig.h"
#include "set.h"


typedef struct _params
{
	char* bam_file; /* paths to comma separated input BAM files as a single string before being tokenized */
	char* outprefix; /* prefix for the output files */
	char *sonic_file;
	char *logfile;
	sv_type svs_to_find;
	unsigned int threads; 
	int chromosome_count;	 
    int barcode_len;    //
    int max_frag_size;  //
    int max_support;
    int sample_size;
    int min_qual;
    int ploidy;
    int *chr_copy_count; 
    double quasi_clique_lambda;
    double quasi_clique_gamma;
    vector_t *haplotype_chrs;   
    int coverage;
    _Bool filter_gap;  //
    _Bool filter_satellite; //
	_Bool low_mem;
} parameters;


/* Parameter related VALOR functions */
parameters *get_params(void);

void free_params(void /*parameters*/ *vp);

parameters *parse_args(int argc, char **argv);

#endif
