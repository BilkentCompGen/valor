#include "clique_inter.h"

#include <limits.h>
#include <assert.h>
#define INITIAL_CLIQUE_COUNT 12


int icclique_plateau_move( icclique_t *c, graph_t *sv_graph, adjlist_t *al, int TABU){
	double max;
	int added, removed;
	double score;
	ic_sv_t *best_sv = NULL;
	size_t best_index = -1;
	max = INT_MIN;
	int i;
	added = 1;
	removed = 1;
	ic_sv_t *tmp;
	ic_sv_t **etmp;

	//ADD
	for(i=0;i< al->size;i++){
		tmp = al_get_value(al,i);
		if(tmp->inactive){continue;}
		if( tmp->tabu > 0){
			tmp->tabu--;
		}
		if( tmp->tabu == 0 && tmp->dv > c->gamma * (c->v_prime - 1) &&
			(c->e_prime + tmp->dv >= (c->lambda * (c->v_prime * (c->v_prime-1)))/2.0)){
			score = (2.0* (c->e_prime + tmp->dv)) / (c->v_prime * (c->v_prime +1));
			if(max < score){
				best_sv = tmp;
				max = score;
				best_index = i;
			}
		}
	}
	if( max == INT_MIN){ added = 0;}
	else{
		icclique_add_node(c,sv_graph, best_sv);

		best_sv->inactive = 1;
		bucket_t *edges = al_get_edges(al,best_index);
		for(i=0;i<edges->size;i++){
			etmp = bucket_get(edges,i);
			(*etmp)->dv++;
		}
		best_sv->tabu = TABU;
	}

	//REMOVE
	max = INT_MIN;
	for(i=0;i< c->items->size;i++){
		tmp = *(ic_sv_t **)vector_get(c->items,i);
		if( tmp->tabu > 0){
			tmp->tabu--;
		}

		if( tmp->tabu == 0 && tmp->dv < c->gamma * (c->v_prime - 1) &&
			(c->e_prime - tmp->dv > (c->lambda * ((c->v_prime-2) * (c->v_prime-1)))/2.0)){
			score = (2.0* (c->e_prime - tmp->dv)) / ((c->v_prime-1) * (c->v_prime -2));
			if(max < score){
				best_sv = tmp;
				max = score;
				best_index = i;
			}
		}
	}

	if(max==INT_MIN){removed = 0;}
	else{
		icclique_rem_node(c,sv_graph, best_sv);
		best_sv->tabu = TABU;
		best_sv->inactive = 0;
		bucket_t *rem_edges = graph_get_edges(sv_graph,best_sv);
		for(i=0;i<rem_edges->size;i++){
			ic_sv_t **_tmp = bucket_get(rem_edges,i);
			(*_tmp)->dv--;
		}
	}

	return added||removed;
}


icclique_t *icclique_find_icclique(graph_t *sv_graph, vector_t *component, int seed, float lambda, float gamma){

	//Sort nodes with degree
	//qsort(sv_graph->nodes->items, sv_graph->nodes->size, sizeof(void *),graph_node_cmp);
	adjlist_t *adjg = graph_make_al(sv_graph,component);
	if(adjg->size == 0){
		return  NULL;
	}
	al_sortbydegree(adjg);
	
	ic_sv_t *initial = al_get_value(adjg,0);	

	//graph_node *initial = vector_head(sv_graph->nodes);
	icclique_t *icclique = icclique_init(lambda,gamma,initial,sv_graph->number_of_items);
	if(icclique==NULL){ 
//		fprintf(stderr,"Null Clique");
		return NULL;
	}


	bucket_t *edges = graph_get_edges(sv_graph,initial);
	if(edges ==NULL){
//		fprintf(stderr,"Null edges");
		return icclique;
	}
	if(edges->size==0){
//		fprintf(stderr,"No edges\n");
		return icclique;
	}

	int i;
	for(i=0;i<edges->size;i++){
		ic_sv_t **temp = bucket_get(edges,i);
		(*temp)->dv = 1;
	}
	initial->inactive = 1;




	while(icclique_plateau_move(icclique,sv_graph,adjg,3+component->size/100));
	vector_free(adjg);
	return icclique;
}


//TODO test methods
icclique_t *icclique_init(float lambda, float gamma, ic_sv_t *node, size_t gsize){
	if( node==NULL){ return NULL;}
	icclique_t *new_icclique = getMem(sizeof(icclique_t));
	new_icclique->items = vector_init(sizeof(ic_sv_t*),INITIAL_CLIQUE_COUNT);
	new_icclique->check_set = set_init( gsize*2, sizeof(ic_sv_t *));
	new_icclique->check_set->key_cmp = &interc_sv_compd;
	new_icclique->e_prime = 0;
	new_icclique->v_prime = 1;
	new_icclique->lambda = lambda;
	new_icclique->gamma = gamma;
	vector_put(new_icclique->items,&node);
	return new_icclique;
}

void icclique_free(icclique_t *c){
	if(c==NULL){return;}
	vector_free(c->items);
	set_free(c->check_set);
	freeMem(c,sizeof(c));
}
int icclique_add_node(icclique_t *c, graph_t *g, ic_sv_t *node){

	if(!graph_have_node(g,node)){return -1;}
	bucket_t *edges = graph_get_edges(g,node);
	vector_put(c->items,&node);
	set_put(c->check_set,&node);
	c->v_prime++;
	int i;
	ic_sv_t **adj_item;
	for(i=0;i<edges->size;i++){
		adj_item = bucket_get(edges,i);
//		if(vector_contains(c->items,adj_item)!=-1){
		if(set_has(c->check_set,adj_item)){
			c->e_prime++;
		}
	}
	return c->e_prime;
}
int icclique_rem_node(icclique_t *c, graph_t *g, ic_sv_t *node){
	int i = 0;
	assert(graph_have_node(g,node));
	assert(vector_contains(c->items,&node)!=-1);
	if(!graph_have_node(g,node)){return -1;}

	size_t index = vector_contains(c->items, &node);
	if(index!=-1){
		vector_remove(c->items,index);
		set_remove(c->check_set,&node);
	}else{
		return -1;
	}
	c->v_prime--;
	bucket_t *edges = graph_get_edges(g,node);
	ic_sv_t **adj_item;
	for(i=0;i<edges->size;i++){
		adj_item = bucket_get(edges,i);
//		if(vector_contains(c->items,adj_item)!=-1){
		if(set_has(c->check_set,adj_item)){
			c->e_prime--;
		}
	}
	return c->e_prime;
}

