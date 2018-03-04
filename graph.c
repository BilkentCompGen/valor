#include "graph.h"

graph_t *graph_init(size_t init_size, size_t node_size){
        graph_t *new_graph = ht_init(init_size, node_size, sizeof(vector_t));
        new_graph->val_rmv = &vector_free;
	return new_graph;
}

void graph_put_node(graph_t *g, void *item){
        vector_t *edges = ht_put(g,item);
	memset(edges,0,sizeof(vector_t));
	edges->item_sizeof = sizeof(void**);
	edges->items = (void **) getMem(sizeof(void *) * INIT_EDGE_COUNT);
	edges->limit = INIT_EDGE_COUNT;
	edges->size = 0;
	edges->rmv = &free;
	edges->fragmental = 0;
	edges->REMOVE_POLICY = 0;
}

void graph_put_edge(graph_t *g, void *i1, void *i2){

        vector_t *e1 = ht_get_value(g,i1);
	void *e2 = ht_get(g,i2)->key;
        if(e1==NULL){
                fprintf(stderr,"NULL NODE e1\n");
        }
//TODO: Fix contains do contains check
        vector_put(e1,&e2);
}

void graph_free(graph_t *g){
/*	int i,j,k;
	vector_t *bucket;
	pair_t *pair;
	vector_t *vect;
	for(i=0;i<g->size;i++){
		bucket = g->buckets[i];
		for(j=0;j<bucket->size;j++){
			pair = vector_get(bucket,j);
			vect = pair->value;
			for(k=0;k<vect->size;k++){
				free(vect->items[k]);
			}
			free(vect->items);
		}	
	}
*/
	ht_free(g);
}

int graph_have_node(graph_t *g, void *item){
	return ht_get(g,item)!=NULL;
}

vector_t *graph_get_edges(graph_t *g, void *item){
	return ht_get_value(g,item);
}

void graph_trim(graph_t *g){
	int i;
	adjlist_t *al = graph_to_al(g);
	for(i=0;i<al->size;i++){
		vector_t *edges = al_get_edges(al,i);
		if(edges->size ==0){
			void *item = al_get_value(al,i);
			graph_remove_node(g,item,G_REMOVE_SOFT);
		}
	}
	vector_free(al);
}

int graph_remove_node(graph_t *g, void *item, int hard){
	if(!graph_have_node(g,item)){return -1;}

	int i,j;
	adjlist_t *adj;	vector_t *neigh;
	vector_t *edges;
	switch(hard){
	case G_REMOVE_HARD:
		adj= graph_to_al(g);
		for(i=0;i<adj->size;i++){
			edges = al_get_edges(adj,i);
			int idx = vector_contains(edges,&item);
			if(idx!=-1){
				vector_remove(edges,idx);
			}
		}
		vector_free(adj);
	break;
	case G_REMOVE_UNDIRECTED:
		neigh = graph_get_edges(g,item);
		for(i=0;i<neigh->size;i++){
			edges = graph_get_edges(g,*(void **)vector_get(neigh,i));
			if(edges==NULL){continue;}
			for(j=0;j<edges->size;j++){
				if(memcmp(item,*(void**)vector_get(edges,j),g->key_size)==0){
//					printf("found edges at %d-%d\n",i,j);
					vector_remove(edges,j);
					break;
				}

			}
			//assert(j!=edges->size);
		}
	break;
	case G_REMOVE_SOFT:
	break;
	}

	ht_remove(g,item);
	return 0;
}

void graph_set_rem_function(graph_t *g, void (*rmv)(void*)){
	g->key_rmv = rmv;
}
int adj_degree_comp(const void *i1, const void *i2){
	const pair_t * const *pp1 = i1;
	const pair_t * const *pp2 = i2;
	const pair_t *p1 = *pp1;
	const pair_t *p2 = *pp2;
	const vector_t *v1 = p1->value;
	const vector_t *v2 = p2->value;
	if( v1->size < v2->size) return 1;
	if( v1->size > v2->size) return -1;
	return 0;
}


//This can be written as a macro
void al_sortbydegree(adjlist_t *g){
	qsort(g->items,g->size,sizeof(void *),adj_degree_comp);
}

void *al_get_value(adjlist_t *g, size_t index){
	pair_t *p = vector_get(g,index);
	return p!=NULL?p->key:NULL;
}

vector_t *al_get_edges(adjlist_t *g, size_t index){
	pair_t *p = vector_get(g,index);
	return p!=NULL?p->value:NULL;
}

graph_iter_t *make_graph_iter(graph_t *g){
	return make_ht_iterator(g);
}

graph_iter_t *graph_iter_copy(graph_iter_t *iter){
	graph_iter_t *new_iter = getMem(sizeof(graph_iter_t));
	memcpy(new_iter,iter,sizeof(graph_iter_t));
	return new_iter;
}
int graph_iter_has_next(graph_iter_t *iter){
	return ht_iter_has_next(iter);
}
int graph_iter_next(graph_iter_t *iter){
	return ht_iter_next(iter);
}
void *graph_iter_get_value(graph_iter_t *iter){
	return ht_iter_get_key(iter);
}
vector_t *graph_iter_get_edges(graph_iter_t *iter){
	return ht_iter_get_value(iter);
}
void graph_iter_free(graph_iter_t *iter){
	freeMem(iter,sizeof(graph_iter_t));
}

