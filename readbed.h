#ifndef __READPCS
#define __READPCS
#include "common.h"
#include <htslib/sam.h>
#include <htslib/hts.h>
//#include "processbam.h"
#include "interval10X.h"
#include "vector.h"
vector_t **read_discordant_bed(char *bed_path);
vector_t **read_barcoded_bed(char *bed_path);
vector_t **read_pcs_bed(char *bed_path, int pool_no);
void destroy_beds( vector_t **reads);
#endif

