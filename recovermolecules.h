#ifndef __RECOV_MOLECULES
#define __RECOV_MOLECULES

#include "interval10X.h"
#include "vector.h"
#include "common.h"
#include "sonic/sonic.h"

vector_t **read_molecules_from_bed(char *filename);
void append_molecules_to_bed(vector_t *, char *filename);
void filter_molecules(vector_t *,sonic *,int chr);
int interval_comp(const void *v1, const void *i2);
vector_t *recover_molecules(vector_t *vector);
#endif
