#ifndef __GRAPH
#define __GRAPH
#include "common.h"
#include "vector.h"
#include "hashtable.h"


#define G_REMOVE_UNDIRECTED 2
#define G_REMOVE_HARD 1
#define G_REMOVE_SOFT 0

#define INIT_EDGE_COUNT 8

typedef hashtable_t graph_t;
typedef vector_t adjlist_t;
typedef ht_iter_t graph_iter_t;


graph_iter_t *make_graph_iter(graph_t *g);
graph_iter_t *graph_iter_copy(graph_iter_t *g);
int graph_iter_has_next(graph_iter_t *iter);
int graph_iter_next(graph_iter_t *iter);
void *graph_iter_get_value(graph_iter_t *iter);
vector_t *graph_iter_get_edges(graph_iter_t *iter);
void graph_iter_free(graph_iter_t *iter);

graph_t *graph_init(size_t init_size, size_t node_size);
void graph_put_node(graph_t *g, void *item);
void graph_put_edge(graph_t *g, void *i1, void *i2);
void graph_free(graph_t *g);
void graph_trim(graph_t *g);

int graph_remove_node(graph_t *g, void *item, int hard);

void graph_print(graph_t *g,FILE *);
int graph_have_node(graph_t *g, void *item);
int graph_have_node_wcmp(graph_t *g, void *item, int (*cmp)(const void *, const void *));
vector_t *graph_get_edges(graph_t *g, void *item);
void graph_set_rem_function(graph_t *g, void (*rmv)(void*));
#define graph_to_al(G) ht_to_vector((G))

#define graph_make_al(G,V) ht_select_pairs((G),(V))
#define graph_make_al_wcmp(G,V,cmp) ht_select_pairs_wcmp((G),(V),(cmp))
void al_sortbydegree(adjlist_t *g);
int adj_degree_comp(const void *,const void *);
void *al_get_value(adjlist_t *g, size_t index);
vector_t *al_get_edges(adjlist_t *g, size_t index);
#endif
