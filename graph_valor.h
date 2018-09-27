#ifndef _VALOR_GRAPH
#define _VALOR_GRAPH
#include "structural_variation.h"
#include "vector.h"
typedef struct sv_node{
	sv_t *sv;
	vector_t *adj;
}sv_node;

typedef vector_t sv_graph_t;

int sv_graph_adj_comp(const void *, const void *);
vector_t *dfs_components(sv_graph_t *);
vector_t *call_tarjan(sv_graph_t *);
void sv_graph_sort_by_degree(sv_graph_t *);
sv_t *sv_graph_get_value(sv_graph_t*,int);
vector_t *sv_graph_get_edges(sv_graph_t*,int);
void sv_node_free(void *);
sv_graph_t *sv_graph_init();
void sv_graph_add_node(sv_graph_t *, sv_t *node);
void sv_graph_add_edge(sv_graph_t *, int i, int j);
void sv_graph_remove_node(sv_graph_t *, int i);
void sv_graph_remove_edge(sv_graph_t *, int i, int j);
void sv_graph_free(sv_graph_t *);
sv_graph_t *make_sv_graph_t(vector_t *svs);
#endif
