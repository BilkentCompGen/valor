#ifndef __COPY_NUMBER_VARIANTS__
#define __COPY_NUMBER_VARIANTS__
#include "common.h"
#include "sonic/sonic.h"
#include "structural_variation.h"
#include <float.h>
#include "valorconfig.h"
#include "readbam.h"
double calculate_chromosome_RD_mean(bam_info *,sonic *, int chr_no);
void calculate_GC_histogram(bam_info *, sonic *,int chr_no);
//duplication_t *check_cnv(bam_info *, sonic*, int chr, inversion_t *pair);

#endif
