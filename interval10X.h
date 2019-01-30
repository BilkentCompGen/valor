#ifndef __INTERVAL_10X__
#define __INTERVAL_10X__
#include "vector.h"
#define I10X_VECTOR_GET(V,I) ((interval_10X *)vector_get((V),(I)))
#define IDIS_VECTOR_GET(V,I) ((interval_discordant *)vector_get((V),(I)))
typedef struct _interval10X{
	int start;
	int end;
	unsigned long barcode;
} interval_10X;

typedef struct _interval_pair{
	int start1;
	int end1;
	int start2;
	int end2;
	unsigned long barcode;
	
} interval_pair;
typedef interval_pair interval_discordant;

typedef interval_pair splitmolecule_t;
int interval_can_pair(interval_10X *i1, interval_10X *i2);

int interval_pair_overlaps(interval_pair *p1,interval_pair *i2, int relaxation);

int interval_size(interval_10X *interval);

int interval_overlaps( interval_10X *interval, interval_10X *that, int relaxation);

int interval_distance( interval_10X *interval, interval_10X *that);

int in_range(int point1, int point2, int range);

interval_10X *interval_intersect(interval_10X* interval, interval_10X *that);

void interval_pair_intersect(interval_pair *interval, interval_pair *that);

unsigned long encode_ten_x_barcode(unsigned char* source);

int interval_start_comp(const void *v1, const void *v2);

int interval_comp(const void *v1, const void *v2);

int barcode_comp(const void *v1, const void *v2);

int interval_pair_comp(const void *v1, const void *v2);



//Distance for inversions
int i_distance(int start1, int start2, int end1, int end2);


//Distance for duplications
int d_distance(int start1, int start2, int end1, int end2);

int interval_outer_distance(interval_10X a, interval_10X b);

int interval_inner_distance(interval_10X a, interval_10X b);

splitmolecule_t *splitmolecule_copy(splitmolecule_t *to_copy);
void splitmolecule_destroy(splitmolecule_t *molecule);

vector_t *discover_split_molecules(vector_t *regions);
splitmolecule_t *splitmolecule_init(interval_10X *i1,interval_10X *i2);

size_t scl_binary_search(vector_t *intervals, splitmolecule_t *key);
#endif
