#ifndef INTER_C_SV
#define INTER_C_SV
#include "vector.h"
#include "readbam.h"
#include "interval10X.h"
#include "common.h"
#include "valorconfig.h"
typedef interval_pair splitmolecule_t;
typedef struct inter_split_molecule{
	int start1;
	int start2;
	int end1;
	int end2;
	int chr1;
	int chr2;	
    unsigned long barcode;
} inter_split_molecule_t;

// inter-chrosomal svs
typedef struct ic_sv{
	inter_split_molecule_t AB;
	inter_split_molecule_t CD;
	splitmolecule_t EF;
	int chr_source :16;
	int chr_target :16;
	sv_type type;
	unsigned char supports[3];
	short tabu;
	int dv :30;
	int covered :1;
	int inactive :1;
}ic_sv_t;



#define INTERC_BACK_COPY 4
#define INTERC_FORW_COPY 3

void ic_sv_bed_print(FILE *stream, ic_sv_t *sv);


size_t barcode_binary_search(vector_t *mols, unsigned long key);

int inter_split_overlaps(inter_split_molecule_t s1, inter_split_molecule_t s2);

inter_split_molecule_t *inter_split_init(barcoded_read_pair *pair, interval_10X *a, interval_10X *b);

int inter_split_indicates_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2);

ic_sv_t *inter_sv_init(inter_split_molecule_t *s1, inter_split_molecule_t *s2, splitmolecule_t *tra_del ,sv_type type, int orient);

vector_t *find_direct_translocations(vector_t *sp1, vector_t *sp2, vector_t **molecules);
int is_inter_chr_split(barcoded_read_pair *pair, interval_10X *a, interval_10X *b);


// returns a vector of inter_split_molecules
vector_t *find_separated_molecules(vector_t *reads, vector_t *mol_a, vector_t *mol_b);



vector_t *find_interchromosomal_events(vector_t **molecules, bam_vector_pack **reads);
#endif
