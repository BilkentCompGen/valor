#ifndef __CLUSTER
#define __CLUSTER
#include "structural_variation.h"
#include "clique.h"
#include "graph.h"
#include "interval10X.h"
typedef struct __sv_cluster{
        vector_t *items;
        interval_pair *break_points;
	int supports[2];
} sv_cluster;

int cluster_comp(const void *a, const void *b);

sv_cluster *sv_cluster_make(clique_t *c);
void sv_cluster_destroy(void *cluster);
void sv_cluster_graph_fix(sv_cluster *cluster, graph_t *graph);//Remove Cluster from the graph
#endif
