#include "cluster.h"

int cluster_comp( const void *a, const void *b){
	sv_cluster *c1 = *(sv_cluster **)a;
	sv_cluster *c2 = *(sv_cluster **)b;
	return - ((c1->supports[1] + c1->supports[0]) * c1->items->size - (c2->supports[1] + c2->supports[0]) * c2->items->size);
}


sv_cluster *sv_cluster_make(clique_t *c){
	if(c->items->size == 0){return NULL;}
	sv_cluster *new_cluster = getMem(sizeof(sv_cluster));
	new_cluster->items = vector_init(sizeof(sv_t),c->items->size);
	new_cluster->supports[0] = 0;
	new_cluster->supports[1] = 0;

	//	inversion_t **inv_tmp;
	sv_t **sv_tmp;
	int i;
	splitmolecule_t *scl_tmp;


	sv_t *isv = *(sv_t **)vector_get(c->items,0); //Initial Inversion
	//inversion_t *iinv = *(inversion_t **)vector_get(c->items,0); //Initial Inversion

	vector_put(new_cluster->items,isv);

	new_cluster->supports[0]+=(isv)->supports[0];
	new_cluster->supports[1]+=(isv)->supports[1];

	new_cluster->break_points = sv_reduce_breakpoints(isv);
	for(i=1;i<c->items->size;i++){

		sv_tmp = vector_get(c->items,i);
		scl_tmp = sv_reduce_breakpoints(*sv_tmp);

		if(!interval_pair_overlaps(new_cluster->break_points,scl_tmp,CLONE_MEAN)){splitmolecule_destroy(scl_tmp); continue;}
		vector_put(new_cluster->items,*sv_tmp);	
		interval_pair_intersect(new_cluster->break_points, scl_tmp);

		new_cluster->supports[0]+=(*sv_tmp)->supports[0];
		new_cluster->supports[1]+=(*sv_tmp)->supports[1];

		splitmolecule_destroy(scl_tmp);
	}
	return new_cluster;
}


void sv_cluster_destroy(void *vcluster){
	sv_cluster *cluster = vcluster;
	vector_free(cluster->items);
	freeMem(cluster->break_points,sizeof(interval_pair));
	freeMem(cluster,sizeof(sv_cluster));
}

void sv_cluster_graph_fix(sv_cluster *c, vector_t *component, graph_t *g){
	
	qsort(component->items,component->size,sizeof(sv_t *), &sv_comp);
	qsort(c->items->items,c->items->size,sizeof(sv_t *), &sv_comp);

	g_remove_all(g,component,c->items);

/*
	int i;


	graph_iter_t *it = make_graph_iter(g);

	while(graph_iter_has_next(it)){
		sv_t *to_add = graph_iter_get_value(it);
		if(to_add->inactive) {
			graph_iter_next(it);
			continue;
		}
		splitmolecule_t *scl = sv_reduce_breakpoints(to_add);
		if(interval_pair_overlaps(scl,c->break_points,CLONE_MEAN)){
			vector_put(c->items,to_add);
			interval_pair_intersect(c->break_points,scl);
			c->supports[1]+=to_add->supports[1];
			c->supports[0]+=to_add->supports[0];
		}
		splitmolecule_destroy(scl);
		graph_iter_next(it);
	}
	graph_iter_free(it);

	
	for(i=0;i<c->items->size;i++){
		sv_t *to_remove = vector_get(c->items,i);
		graph_remove_node(g,to_remove,G_REMOVE_UNDIRECTED);
	}
*/
}
