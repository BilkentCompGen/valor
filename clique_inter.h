#ifndef _CLIQUE_H
#define _CLIQUE_H
#include "vector.h"
#include "set.h"
#include "graph.h"
#include "common.h"
#include "valorconfig.h"
#include "interc_sv.h"
typedef struct __icclique_t{
	vector_t *items; //Vector of Nodes
	set_t *check_set;
	int e_prime;
	int v_prime;
	float lambda;
	float gamma;
} icclique_t;



icclique_t *icclique_find_icclique(graph_t *sv_graph, vector_t* component, int seed ,float lambda, float gamma);
icclique_t *icclique_init( float lambda, float gamma, ic_sv_t *node,size_t csize);
int icclique_add_node(icclique_t *c, graph_t *g, ic_sv_t *node);
int icclique_rem_node(icclique_t *c, graph_t *g, ic_sv_t *node);

void icclique_free(icclique_t *c);
#endif
