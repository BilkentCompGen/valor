#ifndef __COPY_NUMBER_VARIANTS__
#define __COPY_NUMBER_VARIANTS__
#include "common.h"
#include "sonic/sonic.h"
#include "structural_variation.h"
#include <float.h>
#include "valorconfig.h"
#include "readbam.h"

double get_depth_region(short *depth, int start, int end);
double get_depth_deviation(short *depth, int start, int end);
double make_global_molecule_mean(short *depths, sonic *snc,int chr);
double make_global_molecule_std_dev(short *depths, sonic *snc, int chr, double mean);
short *make_molecule_depth_array(vector_t *regions, sonic *, int chr);
//double calculate_chromosome_RD_mean(bam_info *,sonic *, int chr_no);
//void calculate_GC_histogram(bam_info *, sonic *,int chr_no);
//duplication_t *check_cnv(bam_info *, sonic*, int chr, inversion_t *pair);

#endif
