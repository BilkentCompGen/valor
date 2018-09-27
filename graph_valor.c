#include "graph_valor.h"


int  sv_graph_adj_comp(const void *v1, const void *v2){

	const vector_t *ve1 = (*(sv_node **)v1)->adj;
	const vector_t *ve2 = (*(sv_node **)v2)->adj;


	if(ve1->size < ve2->size){return 1;}
	if(ve1->size > ve2->size){return -1;}
	return 0;
}

void dfs_step(sv_graph_t *g, vector_t *comp, int u, unsigned char* visited){
	visited[u] = 1;
	vector_put(comp,&u);
	int i;
	vector_t *edges = sv_graph_get_edges(g,u);
	for(i=0;i<edges->size;i++){
		if(visited[i]!=1){
			dfs_step(g,comp,i,visited);
		}
	}
}

vector_t *dfs_components(sv_graph_t *g){
	vector_t *comps = vector_init(sizeof(vector_t),48);
	unsigned char* visited = malloc(sizeof(unsigned char) *g->size);
	int i;
	for(i=0;i<g->size;i++){
		visited[i] = 0;
	}
	for(i=0;i<g->size;i++){
		vector_t *comp = vector_init(sizeof(int),212);
		if(visited[i]==0){
			dfs_step(g,comp,i,visited);
		}
		vector_soft_put(comps,comp);
	}
	free(visited);
	return comps;
}

void tarjan_step(sv_graph_t *g, vector_t *comps,  int u, int *disc, int *low, vector_t *stack, unsigned char* stack_member){
	static int time = 0;
	disc[u] = low[u] = ++time;
	vector_put(stack, &u);
	stack_member[u] = 1;

	int i;
	vector_t *edges= sv_graph_get_edges(g,u);
	for(i=0;i<edges->size;i++){
		if( disc[i] == -1){
			tarjan_step(g,comps,i,disc,low,stack,stack_member);
			
			low[u] = MIN( low[u], low[i]);
		}
		else if (stack_member[i]){
			low[u] = MIN(low[u], disc[i]);
		}

	}
	int w = 0;
	
	if(low[u] == disc[u]){
		int *top = vector_tail(stack);
		vector_soft_put(comps,vector_init(sizeof(int),8));
		while(*top !=u){
			w = *(int *) vector_tail(stack);
			vector_put(vector_tail(comps),&w);
			stack_member[w] = 0;	
			vector_remove(stack,stack->size-1);
			top = vector_tail(stack);
		}
		w = *(int *) vector_tail(stack);
		vector_put(vector_tail(comps),&w);
		stack_member[w] = 0;	
		vector_remove(stack,stack->size-1);
	}
}

vector_t *call_tarjan(sv_graph_t *g){
	int *disc = malloc(sizeof(int) * g->size);
	int *low = malloc(sizeof(int) *g->size);
	unsigned char *stack_member = malloc(sizeof(unsigned char) * g->size);
	int i;
	for(i=0;i<g->size;i++){
		disc[i] = low[i] = -1;
		 stack_member[i] = 0;
	}
	vector_t *stack = vector_init(sizeof(int),1+g->size/2);
	vector_t *comps = vector_init(sizeof(vector_t),16);
	for(i=0;i<g->size;i++){
		if(disc[i]==-1){
			tarjan_step(g,comps,i,disc,low,stack,stack_member);
		}
	}
	vector_free(stack);
	free(disc);
	free(low);
	free(stack_member);
	return comps;
}
int sv_node_adj_cmp( const void *a1, const void *a2){
	const sv_node * const *sv1 = a1;
	const sv_node * const *sv2 = a2;
	const sv_node *s1 = *sv1;
	const sv_node *s2 = *sv2;
	if( s1->adj->size < s2->adj->size) return 1;
	if( s1->adj->size > s2->adj->size) return -1;
	return 0;
}

void sv_graph_sort_by_degree(sv_graph_t *g){
	qsort(g->items,g->size,sizeof(void*),sv_node_adj_cmp);
}	
void sv_node_free(void *vnode){
	sv_node *node = vnode;
	vector_free(node->adj);
	sv_destroy(node->sv);
	free(vnode);
}


sv_t *sv_graph_get_value(sv_graph_t *g, int i){
	return ((sv_node*)vector_get(g,i))->sv;
}

vector_t *sv_graph_get_edges(sv_graph_t *g, int i){
	return ((sv_node*)vector_get(g,i))->adj;
}
sv_graph_t *sv_graph_init(int init_size){
	sv_graph_t *g = vector_init(sizeof(sv_node),init_size);
	g->rmv = &sv_node_free;
	return g;
}

void sv_graph_put_node(sv_graph_t *g,sv_t *node){
	vector_put(g,&(sv_node){.sv=node,.adj=vector_init(sizeof(sv_node *),16)});
	sv_node *tail = vector_tail(g);
	tail->adj->REMOVE_POLICY = REMP_FAST;
}

void sv_graph_put_edge(sv_graph_t *g, int i, int j){
	sv_node *n1 = vector_get(g,i);
	sv_node *n2 = vector_get(g,j);
	vector_put(n1->adj,&n2);
	vector_put(n2->adj,&n1);
}

void sv_graph_remove_node(sv_graph_t *g, int i){

	vector_t *edges= sv_graph_get_edges(g,i);
	
	while(edges->size > 0){
		sv_node **node = vector_get(edges,edges->size-1);
		int jindx = vector_contains(g,*node);
		sv_graph_remove_edge(g,i,jindx);

	}
	vector_remove(g,i);
}

void sv_graph_remove_edge(sv_graph_t *g, int i, int j){
	sv_node *n1 = vector_get(g,i);
	sv_node *n2 = vector_get(g,j);
	int i1 = vector_contains(n1->adj,&n2);
	int i2 = vector_contains(n2->adj,&n1);
	vector_remove(n1->adj,i1);
	vector_remove(n2->adj,i2);
}

void sv_graph_free(sv_graph_t *g){
	vector_free(g);
}


sv_graph_t *make_sv_graph_t(vector_t *svs){
        sv_graph_t *g = sv_graph_init(svs->size);

        int i,j;
        printf("Adding Nodes\n");
        for(i=0;i<svs->size;i++){
                sv_graph_put_node(g,vector_get(svs,i));
        }
        for(i=0;i<svs->size;i++){
                sv_t *a = vector_get(svs,i);
                for(j=i+1;j<svs->size;j++){
                        sv_t *b = vector_get(svs,j);

                        if(sv_overlaps(a,b)){
                                sv_graph_put_edge(g,i,j);

                        }
                }
        }
        return g;
}

