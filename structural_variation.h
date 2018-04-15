#ifndef __INVERSION
#define __INVERSION
#include "vector.h"
#include "graph.h"
#include "interval10X.h"
#include "readbam.h"


typedef interval_pair splitmolecule_t;

typedef struct sv{
	splitmolecule_t AB;
	splitmolecule_t CD;
//	double depths[2];
	int supports[2];
	int tabu;
	int dv;
	sv_type type;
	char orientation;
	_Bool covered;
	_Bool inactive;
}sv_t;

#define DUP_BACK_COPY 2
#define DUP_FORW_COPY 1
#define SPCL_VECTOR_GET(V,I) ((splitmolecule_t *)vector_get((V),(I)))

#define SV_VECTOR_GET(V,I) ((sv_t *)vector_get((V),(I)))

int g_remove_all(graph_t *,vector_t *, vector_t *);
vector_t *g_dfs_components(graph_t *);
const char *sv_type_name(sv_type);

int _svcmp(const void *v1, const void *v2,size_t);
int sv_comp(const void *v1, const void *v2);
int sv_equals(const void *i1, const void *i2);
splitmolecule_t *splitmolecule_init(interval_10X *i1,interval_10X *i2);
void sv_fprint(FILE *stream, int chr, sv_t *inv);
graph_t *make_sv_graph(vector_t *svs);
sv_t *sv_init(splitmolecule_t *,splitmolecule_t *,sv_type);


void sv_destroy(void/*sv_t*/ *dup);
splitmolecule_t *sv_reduce_breakpoints(sv_t*);
vector_t *discover_split_molecules(vector_t *regions);
vector_t *find_svs(vector_t *split_molecules,sv_type);

void sv_graph_reset(graph_t *g);
void splitmolecule_destroy(splitmolecule_t *molecule);
void update_sv_supports(vector_t *svs, bam_vector_pack *,sv_type);
size_t sv_hf(hashtable_t *table, const void *inv);
splitmolecule_t *splitmolecule_copy(splitmolecule_t *to_copy);
void sv_graph_visualize(graph_t *g, FILE *stream);

int inversion_overlaps(sv_t *i1, sv_t *i2);
int duplication_overlaps(sv_t *i1, sv_t *i2);
int sv_overlaps(sv_t *i1, sv_t *i2);
#endif
