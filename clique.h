#ifndef _CLIQUE_H
#define _CLIQUE_H
#include "vector.h"
#include "graph.h"
#include "common.h"
#include "valorconfig.h"
#include "structural_variation.h"
typedef struct __clique_t{
	vector_t *items; //Vector of Nodes
	int e_prime;
	int v_prime;
	float lambda;
	float gamma;
} clique_t;


int clique_test(int argc, char **argv);
clique_t *clique_find_clique(graph_t *sv_graph, int seed ,float lambda, float gamma);
clique_t *clique_init( float lambda, float gamma, sv_t *node);
int clique_add_node(clique_t *c, graph_t *g, sv_t *node);
int clique_rem_node(clique_t *c, graph_t *g, sv_t *node);

void clique_free(clique_t *c);
#endif
