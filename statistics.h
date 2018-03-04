#ifndef __STATISTICS__ 
#define __STATISTICS__ 
#include "common.h"
#include "vector.h"
#include "interval10X.h"
#include "valorconfig.h"
#include <math.h>
double molecule_mean(vector_t *);
double molecule_std(vector_t *, double mean);
vector_t *group_overlapping_molecules(vector_t *);

double molecule_group_std(vector_t *, double mean);
#endif                 

