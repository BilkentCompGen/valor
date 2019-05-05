#include "valorconfig.h"
#include "structural_variation.h"
#include <stdio.h>
#include "progress.h"
#include "sonic/sonic.h"
#include "cnv.h"




void sv_fprint(FILE *stream, int chr, sv_t *t){
	fprintf(stream,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
			chr,
			t->AB.start1,
			t->AB.end1,
			t->AB.start2,
			t->AB.end2,
			t->CD.start1,
			t->CD.end1,
			t->CD.start2,
			t->CD.end2,
			t->supports[0],
			t->supports[1],
			sv_type_name(t->type));
}


int _svcmp(const void *v1, const void *v2, size_t size){
	const sv_t *s1 = v1;
	const sv_t *s2 = v2;
	if(s1->AB.start1 != s2->AB.start2){
		return s1->AB.start1 - s2->AB.start1;
	}
	if(s1->AB.end1 != s2->AB.end1){
		return s1->AB.end1 - s2->AB.end1;
	}
	if(s1->AB.start2 != s2->AB.start2){
		return s1->AB.start2 - s2->AB.start2;
	}
	if(s1->AB.end2 != s2->AB.end2){
		return s1->AB.end2 - s2->AB.end2;
	}
	if(s1->CD.start1 != s2->CD.start2){
		return s1->CD.start1 - s2->CD.start1;
	}
	if(s1->CD.end1 != s2->CD.end1){
		return s1->CD.end1 - s2->CD.end1;
	}
	if(s1->CD.start2 != s2->CD.start2){
		return s1->CD.start2 - s2->CD.start2;
	}
	if(s1->CD.end2 != s2->CD.end2){
		return s1->CD.end2 - s2->CD.end2;
	}
	return 0;
}

int sv_compd(const void *v1, const void *v2, size_t val){
	return sv_comp(v1,v2);
}
int sv_comp(const void *v1, const void *v2){
	sv_t *s1 = *(void **)v1;
	sv_t *s2 = *(void **)v2;
	if(s1->AB.start1 != s2->AB.start2){
		return s1->AB.start1 - s2->AB.start1;
	}
	if(s1->AB.end1 != s2->AB.end1){
		return s1->AB.end1 - s2->AB.end1;
	}
	if(s1->AB.start2 != s2->AB.start2){
		return s1->AB.start2 - s2->AB.start2;
	}
	if(s1->AB.end2 != s2->AB.end2){
		return s1->AB.end2 - s2->AB.end2;
	}
	if(s1->CD.start1 != s2->CD.start2){
		return s1->CD.start1 - s2->CD.start1;
	}
	if(s1->CD.end1 != s2->CD.end1){
		return s1->CD.end1 - s2->CD.end1;
	}
	if(s1->CD.start2 != s2->CD.start2){
		return s1->CD.start2 - s2->CD.start2;
	}
	if(s1->CD.end2 != s2->CD.end2){
		return s1->CD.end2 - s2->CD.end2;
	}
	return 0;
}

//TODO consider barcode
int sv_equals(const void *i1, const void *i2){
	const sv_t *ii1 = i1;
	const sv_t *ii2 = i2;
	return !(ii1->AB.start1 == ii2->AB.start1 &&
			ii1->AB.start2 == ii2->AB.start2 &&
			ii1->AB.end1 == ii2->AB.end1 &&
			ii1->AB.end2 == ii2->AB.end2 &&
			ii1->CD.start1 == ii2->CD.start1 &&
			ii1->CD.start2 == ii2->CD.start2 &&
			ii1->CD.end1 == ii2->CD.end1 &&
			ii1->CD.end2 == ii2->CD.end2 &&
			ii1->AB.barcode == ii2->AB.barcode &&
			ii1->CD.barcode == ii2->CD.barcode &&
			ii1->type == ii2->type);

	//	return 	memcmp(&(ii1->AB),&(ii2->AB),sizeof(splitmolecule_t)) &&
	//		memcmp(&(ii1->CD),&(ii2->CD),sizeof(splitmolecule_t));	
}


void sv_g_dfs_step(graph_t *g, vector_t *comp, sv_t *sv){
	sv->covered = 1;
	vector_put(comp,sv);
	vector_t *edges = graph_get_edges(g,sv);
	int i;
	for(i=0;i<edges->size;i++){
		sv_t **val = vector_get(edges,i);
		if(!(*val)->covered){
			sv_g_dfs_step(g,comp,*val);
		}
	}
}

vector_t *sv_g_dfs_components(graph_t *g){
	adjlist_t *al = graph_to_al(g);
	vector_t *comps = vector_init(sizeof(vector_t),16);
	comps->rmv  = &vector_free;
	int i;

	for(i=0;i<al->size;i++){
		sv_t *sv = al_get_value(al,i);
		if( !sv->covered){
			vector_t *comp = vector_init(sizeof(sv_t),40);
			comp->rmv = &sv_destroy;
			sv_g_dfs_step(g,comp,sv);
			vector_soft_put(comps,comp);
		}

	}
	vector_free(al);
	return comps;
}

#define BIG_PRIME 1300501
size_t sv_hf(hashtable_t *table, const void *vsv){
	const sv_t *sv = vsv;
	size_t hash = 0;
	hash+=sv->AB.barcode;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->CD.barcode;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->AB.start1;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->AB.start2;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->AB.end1;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->AB.end2;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->CD.start1;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->CD.start2;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->CD.end1;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->CD.end2;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	hash+=sv->type;
	hash*=BIG_PRIME;
	hash= hash % table->size;
	return hash;


	//return SuperFastHash((vsv),2*sizeof(splitmolecule_t)) % table->size;
}

int splitmolecule_indicates_inverted_duplication(splitmolecule_t s1, splitmolecule_t s2){

	if( !(s1.start1 < s2.start1 &&
				s1.end1 < s2.end1 &&
				s1.start2 > s2.start2 &&
				s1.end2 > s2.end2)){
		return 0;
	}

	interval_10X A = (interval_10X){s1.start1,s1.end1,s1.barcode};  
	interval_10X B = (interval_10X){s1.start2,s1.end2,s1.barcode};  
	interval_10X C = (interval_10X){s2.start1,s2.end1,s2.barcode};  
	interval_10X D = (interval_10X){s2.start2,s2.end2,s2.barcode};  
	if(interval_inner_distance(A,C) < DUP_GAP &&
			interval_inner_distance(A,C) > DUP_OVERLAP &&
			interval_outer_distance(B,D) < DUP_MAX_SIZE &&
			interval_outer_distance(B,D) > DUP_MIN_SIZE){
		return DUP_BACK_COPY;
	}
	else if(interval_inner_distance(B,D) < DUP_GAP &&
			interval_inner_distance(B,D) > DUP_OVERLAP &&
			interval_outer_distance(A,C) < DUP_MAX_SIZE &&
			interval_outer_distance(A,C) > DUP_MIN_SIZE){
		return DUP_FORW_COPY;
	}
	return 0;
}
//TODO Sanity Check
int splitmolecule_indicates_duplication(splitmolecule_t s1, splitmolecule_t s2){

	if( !(s1.start1 < s2.start1 &&
				s1.end1 < s2.end1 &&
				s1.start2 < s2.start2 &&
				s1.end2 < s2.end2)){
		return 0;
	}

	interval_10X A = (interval_10X){s1.start1,s1.end1,s1.barcode};  
	interval_10X B = (interval_10X){s1.start2,s1.end2,s1.barcode};  
	interval_10X C = (interval_10X){s2.start1,s2.end1,s2.barcode};  
	interval_10X D = (interval_10X){s2.start2,s2.end2,s2.barcode};  
	if(interval_inner_distance(A,C) < DUP_GAP &&
			interval_inner_distance(A,C) > DUP_OVERLAP &&
			interval_outer_distance(B,D) < DUP_MAX_SIZE &&
			interval_outer_distance(B,D) > DUP_MIN_SIZE){
		return DUP_BACK_COPY;
	}
	else if(interval_inner_distance(B,D) < DUP_GAP &&
			interval_inner_distance(B,D) > DUP_OVERLAP &&
			interval_outer_distance(A,C) < DUP_MAX_SIZE &&
			interval_outer_distance(A,C) > DUP_MIN_SIZE){
		return DUP_FORW_COPY;
	}
	return 0;
}

//TODO make sure different barcode
int splitmolecule_indicates_inversion(splitmolecule_t s1, splitmolecule_t s2){
	if( abs(s2.end2-s1.start1)> INV_MAX_SIZE){ return 0;}
	return s1.start1 < s2.start1 &&
		s1.end1 < s2.end1 &&
		s1.start2 < s2.start2 &&
		s1.end2 < s2.end2 &&
		i_distance(s2.start1,s1.start1,s2.end1,s1.end1) > INV_OVERLAP &&
		i_distance(s2.start2,s1.start2,s2.end2,s1.end2) > INV_OVERLAP &&
		i_distance(s2.start1,s1.start1,s2.end1,s1.end1) < INV_GAP &&
		i_distance(s2.start2,s1.start2,s2.end2,s1.end2) < INV_GAP;
}


int splitmolecule_indicates_sv(splitmolecule_t *s1, splitmolecule_t *s2, sv_type type){
	switch(type){
		case SV_INVERSION:
			return splitmolecule_indicates_inversion(*s1,*s2);
		case SV_DIRECT_DUPLICATION:
			return splitmolecule_indicates_duplication(*s1,*s2);
		case SV_INVERTED_DUPLICATION:
			return splitmolecule_indicates_inverted_duplication(*s1,*s2);
        case SV_TRANSLOCATION:
			return splitmolecule_indicates_duplication(*s1,*s2);
        case SV_INVERTED_TRANSLOCATION:
			return splitmolecule_indicates_inverted_duplication(*s1,*s2);
        default:
			fprintf(stderr,"Unknown SV type ordinal %d\n",type);
			VALOR_LOG("Unknown SV type ordinal %d\n",type);
			exit(-1);
	}
	return 0;
}

void sv_reset(sv_t *sv){
	sv->covered = 0;
	sv->tabu = 0;
	sv->dv = 0;
	sv->inactive = 0;
}



sv_t *sv_init(splitmolecule_t *sc1,splitmolecule_t *sc2,sv_type type){
	sv_t *new_i = getMem(sizeof(sv_t));
	memset(new_i,0,sizeof(sv_t));
	new_i->supports[0] = 1;//TODO fix this
	new_i->supports[1] = 1;//TODO fix this

	new_i->covered = 0;
	new_i->tabu = 0;
	new_i->dv = 0;
	new_i->inactive = 0;
	new_i->AB=*sc1;

	if(sc2!=NULL){//for sv's with 1 split molecule
		new_i->CD=*sc2;
	}
	else{
		new_i->CD=(splitmolecule_t){0,0,0,0,0L};
	}
	new_i->type=type;
	return new_i;
}

void sv_destroy(void *sv){
	freeMem(sv,sizeof(sv_t));
}






int inversion_overlaps(sv_t *i1, sv_t *i2){
	return interval_pair_overlaps(
			&(splitmolecule_t){
			MIN(i1->AB.end1,i1->CD.start1),MAX(i1->AB.end1,i1->CD.start1),
			MIN(i1->AB.end2,i1->CD.start2),MAX(i1->AB.end2,i1->CD.start2)
			},
			&(splitmolecule_t){
			MIN(i2->AB.end1,i2->CD.start1),MAX(i2->AB.end1,i2->CD.start1),
			MIN(i2->AB.end2,i2->CD.start2),MAX(i2->AB.end2,i2->CD.start2)
			},CLONE_MEAN);
}

//TODO change this if it needs change
int duplication_overlaps(sv_t *i1, sv_t *i2){
	return interval_pair_overlaps(
			&(splitmolecule_t){
			MIN(i1->AB.end1,i1->CD.start1),MAX(i1->AB.end1,i1->CD.start1),
			MIN(i1->AB.end2,i1->CD.start2),MAX(i1->AB.end2,i1->CD.start2)
			},
			&(splitmolecule_t){
			MIN(i2->AB.end1,i2->CD.start1),MAX(i2->AB.end1,i2->CD.start1),
			MIN(i2->AB.end2,i2->CD.start2),MAX(i2->AB.end2,i2->CD.start2)
			},CLONE_MEAN);
}

int tandem_duplication_overlaps(sv_t *i1, sv_t *i2){
	return interval_pair_overlaps(
			&(i1->AB),&(i2->AB)
			,CLONE_MEAN);
}


int deletion_overlaps(sv_t *i1, sv_t *i2){
	return interval_pair_overlaps(
			&(i1->AB),&(i2->AB)
			,CLONE_MEAN);
}

int sv_overlaps(sv_t *i1, sv_t *i2){
	if(i1->type!=i2->type){ return 0;}
	switch(i1->type){
		case SV_INVERSION:
			return inversion_overlaps(i1,i2);

		case SV_DIRECT_DUPLICATION:
			return duplication_overlaps(i1,i2);
		case SV_INVERTED_DUPLICATION:
			return duplication_overlaps(i1,i2);
		case SV_TANDEM_DUPLICATION:
			return tandem_duplication_overlaps(i1,i2);
        case SV_TRANSLOCATION:
			return duplication_overlaps(i1,i2);
		case SV_INVERTED_TRANSLOCATION:
			return duplication_overlaps(i1,i2);
		case SV_DELETION:
			return deletion_overlaps(i1,i2);

		default:
			fprintf(stderr,"Unknown SV type ordinal %d\n",i1->type);
			VALOR_LOG("Unknown SV type ordinal %d\n",i1->type);
			exit(-1);
	}
	return 0;
}


#define MAGIC_NODE_MARKER -126
/*
 *
 *  This removes every element in the `items` from `g`
 *  Nodes in the `items` shouldn't have any edges outside of the `component`
 *
 *  Achtung: 
 *  	sv->dv shouldn't be -126. Make sure it is not -126 outside of this function
 *  	Make sure `g`->hf does not use sv->dv
 *  	Make sure `g`->key_cmp does not use sv->dv 
 *  		(default is memcmp(a,b,sizeof(sv_t)), so Make sure it is changed)
 */
int g_remove_all(graph_t *g, vector_t *component,vector_t *items){
	int i,j;

	for(i=0;i<items->size;i++){
		void *item = vector_get(items,i);
		//                graph_remove_node(g,item,G_REMOVE_SOFT);
		sv_t *sv = ht_get(g,item)->key;
		sv->dv = MAGIC_NODE_MARKER; //Hax
	}

	for(i=0;i<component->size;i++){
		vector_t *edges = graph_get_edges(g,vector_get(component,i));
		if(edges==NULL){continue;}
		edges->REMOVE_POLICY = REMP_LAZY;
		for(j=0;j<edges->size;j++){
			sv_t **ptr = vector_get(edges,j);

			if((*ptr)->dv==MAGIC_NODE_MARKER){
				vector_remove(edges,j);
			}
		}
		vector_defragment(edges);
		edges->REMOVE_POLICY = REMP_FAST;
	}

	for(i=0;i<items->size;i++){
		void *item = vector_get(items,i);
		graph_remove_node(g,item,G_REMOVE_SOFT);
	}


	component->REMOVE_POLICY = REMP_LAZY;
	for(i=0;i<component->size;i++){
		void *ptr = vector_get(component,i);
		if(!graph_have_node(g,ptr)){
			vector_remove(component,i);
		}
	}
	vector_defragment(component);
	component->REMOVE_POLICY = REMP_SORTED;

	return 0;
}


void sv_graph_reset(graph_t *g){
	int i,j;

	for(i=0;i<g->size;i++){
		for(j=0;j<g->buckets[i]->size;j++){
			pair_t *pair = vector_get(g->buckets[i],j);
			sv_reset(pair->key);
		}
	}
}

graph_t *make_sv_graph(vector_t *svs){
	graph_t *g = graph_init(svs->size * 2, sizeof(sv_t));
	g->hf = &sv_hf;
	g->key_cmp = &_svcmp;
	int i,j;
	for(i=0;i<svs->size;i++){
		graph_put_node(g,vector_get(svs,i));
	}
	for(i=0;i<svs->size;i++){
		sv_t *a = vector_get(svs,i);
		for(j=i+1;j<svs->size;j++){
			sv_t *b = vector_get(svs,j);

			if(sv_overlaps(a,b)){
				graph_put_edge(g,a,b);
				graph_put_edge(g,b,a);
			}
		}
	}
	return g;
}


#define SV_INIT_LIMIT 10000

vector_t *find_svs(vector_t *split_molecules, sv_type type, int chr){
	int i,j,k;
	vector_t *svs;


	int num_threads = 1;

#ifdef _OPENMP
	num_threads =  omp_get_num_threads();
#endif
	if(num_threads == 1){

		svs = vector_init(sizeof(sv_t),SV_INIT_LIMIT);
		splitmolecule_t *AB;
		splitmolecule_t *CD;

		for(i=0;i<split_molecules->size;i++){
			AB = vector_get(split_molecules,i);

			if(type&SV_DELETION){
				sv_t *tmp = sv_init(AB,NULL,SV_DELETION);
			    tmp->chr = chr;
                vector_soft_put(svs,tmp);
			}
			if(type&SV_TANDEM_DUPLICATION){
				sv_t *tmp = sv_init(AB,NULL,SV_TANDEM_DUPLICATION);

			            tmp->chr = chr;

				vector_soft_put(svs,tmp);
			}
			if(!(type & ( SV_INVERSION | SV_DIRECT_DUPLICATION | SV_INVERTED_DUPLICATION | SV_TRANSLOCATION | SV_INVERTED_TRANSLOCATION))){ continue;}
			for(j=0;j<split_molecules->size;j++){
				if(i==j){continue;}
				CD = vector_get(split_molecules,j);
				for(k=SV_INVERSION;k<SV_MAX_ID;k=k<<1){
					if((k&type)==0){ continue;}
					char orient = splitmolecule_indicates_sv(AB,CD,k);
					if(orient){
						sv_t *tmp = sv_init(AB,CD,k);
						tmp->orientation = orient;
			            tmp->chr = chr;
                        vector_soft_put(svs,tmp);
					}
				}
			}
		}
	}
	else{

		vector_t **osvs = malloc(sizeof(vector_t *) * num_threads);
		for(i=0;i < num_threads;i++){
			osvs[i]=vector_init(sizeof(sv_t),SV_INIT_LIMIT);
		}

#pragma omp parallel for
		for(i=0;i<split_molecules->size;i++){
			splitmolecule_t *AB = vector_get(split_molecules,i);
			if(type & SV_DELETION){
				sv_t *tmp = sv_init(AB,NULL,SV_DELETION);
		
			            tmp->chr = chr;
                vector_soft_put(osvs[omp_get_thread_num()],tmp);
			}
            if(type&SV_TANDEM_DUPLICATION){
				sv_t *tmp = sv_init(AB,NULL,SV_TANDEM_DUPLICATION);

			            tmp->chr = chr;
                vector_soft_put(osvs[omp_get_thread_num()],tmp);
			}

			if(!(type & ( SV_INVERSION | SV_DIRECT_DUPLICATION | SV_INVERTED_DUPLICATION| SV_TRANSLOCATION | SV_INVERTED_TRANSLOCATION))){ continue;}

			for(j=0;j<split_molecules->size;j++){
				if(i==j){continue;}
				splitmolecule_t *CD = vector_get(split_molecules,j);
				for(k=SV_INVERSION;k<SV_MAX_ID;k=k<<1){
					if((k&type)==0){ continue;}
					char orient = splitmolecule_indicates_sv(AB,CD,type);
					if(orient){
						sv_t *tmp = sv_init(AB,CD,type);
						tmp->orientation = orient;
						
			            tmp->chr = chr;
                        vector_soft_put(osvs[omp_get_thread_num()],tmp);
					}
				}
			}
		}
		size_t vsize = 0;
		for(i=0;i<num_threads;i++){
			vsize+=osvs[i]->size;
		}
		svs = vector_init(sizeof(sv_t),vsize);
		for(i=0;i<num_threads;i++){
			for(j=0;j<osvs[i]->size;j++){
				vector_soft_put(svs,vector_get(osvs[i],j));
			}
			vector_tabularasa(osvs[i]);
			vector_free(osvs[i]);
		}
		free(osvs);
	}

	return svs;
}


size_t scl_binary_searchdup(vector_t *intervals, splitmolecule_t *key){
	if(intervals->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = intervals->size - 1;
	int counter = 0;	
	while( first < last){
		mid = (first + last)/2;
		if(IDIS_VECTOR_GET(intervals,mid)->end1 < key->start1){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
		counter ++;
	}
	return mid;
}




void update_tandem_duplication_supports_b(sv_t *dup, vector_t *mp_reads){
	int j;
	int mp_support = 0;
	int midAB = scl_binary_search(mp_reads,&(dup->AB));

	for(j=midAB;j<mp_reads->size;j++){
		if(interval_pair_overlaps(&(dup->AB),vector_get(mp_reads,j),MAX_FRAG_SIZE)){
			mp_support++;				
		}
		if(dup->AB.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
			break;
		}
		if(mp_support > MAX_SUPPORT){ break;}	
    }
	dup->supports[1] = mp_support;
	dup->supports[0] = mp_support;
}



void update_deletion_supports_b(sv_t *del, vector_t *pm_reads){
	int j;
	int pm_support;
	int midAB;

	pm_support = 0;
	midAB = scl_binary_search(pm_reads,&(del->AB));

	for(j=midAB;j<pm_reads->size;j++){
		if(interval_pair_overlaps(&(del->AB),vector_get(pm_reads,j),MAX_FRAG_SIZE)){//CLONE_MEAN)){
			pm_support++;				
		}
		if(del->AB.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
			break;
		}
		if(pm_support > MAX_SUPPORT){ break;}	
	}



	del->supports[1] = pm_support;
	del->supports[0] = pm_support;
}


void update_deletion_supports(vector_t *dels, vector_t *pm_reads){
	int i,j;
	int pm_support;
	int midAB;

	qsort(pm_reads->items,pm_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	for(i=0;i<dels->size;i++){
		pm_support = 0;
		midAB = scl_binary_search(pm_reads,&(SV_VECTOR_GET(dels,i)->AB));

		for(j=midAB;j<pm_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(dels,i)->AB),vector_get(pm_reads,j),CLONE_MEAN)){
				pm_support++;				
			}
			if(SV_VECTOR_GET(dels,i)->AB.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
				break;
			}
			if(pm_support > MAX_SUPPORT){ break;}	
		}



		SV_VECTOR_GET(dels,i)->supports[1] = pm_support;
		SV_VECTOR_GET(dels,i)->supports[0] = pm_support;
	}

}


void update_duplication_supports_b(sv_t *dup, vector_t *pm_reads, vector_t *mp_reads){
	int j;
	int pm_support;
	int mp_support;
	int midAB;
	int midCD;


	pm_support = 0;
	mp_support = 0;


	if(dup->orientation == DUP_FORW_COPY){
		midAB = scl_binary_search(pm_reads,&(dup->CD));
		midCD = scl_binary_search(mp_reads,&(dup->AB));
		for(j=midAB;j<pm_reads->size;j++){
			//			printf("midAB: %d, j %d",midAB,j);
			if(interval_pair_overlaps(&(dup->CD),vector_get(pm_reads,j),CLONE_MEAN)){
				pm_support++;				
			}
			//			printf("\n");
			if(dup->CD.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
				break;
			}
			if(pm_support > MAX_SUPPORT){ break;}
		}

		for(j=midCD;j<mp_reads->size;j++){
			if(interval_pair_overlaps(&(dup->AB),vector_get(mp_reads,j),CLONE_MEAN)){
				mp_support++;				
			}
			if(dup->AB.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
				break;
			}
			if(mp_support > MAX_SUPPORT){ break;}
		}
	}
	else if(dup->orientation == DUP_BACK_COPY){

		midAB = scl_binary_search(pm_reads,&(dup->AB));
		midCD = scl_binary_search(mp_reads,&(dup->CD));
		for(j=midAB;j<pm_reads->size;j++){
			if(interval_pair_overlaps(&(dup->AB),vector_get(pm_reads,j),CLONE_MEAN)){
				pm_support++;				
			}
			if(dup->AB.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
				break;
			}
			if(pm_support > MAX_SUPPORT){ break;}	
		}

		for(j=midCD;j<mp_reads->size;j++){
			if(interval_pair_overlaps(&(dup->CD),vector_get(mp_reads,j),CLONE_MEAN)){
				mp_support++;				
			}
			if(dup->CD.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
				break;
			}

			if(mp_support > MAX_SUPPORT){ break;}	
		}
	}
	dup->supports[0] = mp_support;
	dup->supports[1] = pm_support;
}



void update_duplication_supports(vector_t *dups, vector_t *pm_reads, vector_t *mp_reads){
	int i,j;
	int pm_support;
	int mp_support;
	int midAB;
	int midCD;

	qsort(pm_reads->items,pm_reads->size,sizeof(interval_discordant *),interval_pair_comp);
	qsort(mp_reads->items,mp_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	for(i=0;i<dups->size;i++){
		pm_support = 0;
		mp_support = 0;


		if(SV_VECTOR_GET(dups,i)->orientation == DUP_FORW_COPY){
			midAB = scl_binary_search(pm_reads,&(SV_VECTOR_GET(dups,i)->CD));
			midCD = scl_binary_search(mp_reads,&(SV_VECTOR_GET(dups,i)->AB));
			for(j=midAB;j<pm_reads->size;j++){
				//			printf("midAB: %d, j %d",midAB,j);
				if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->CD),vector_get(pm_reads,j),MAX_FRAG_SIZE)){
					pm_support++;				
				}
				//			printf("\n");
				if(SV_VECTOR_GET(dups,i)->CD.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
					break;
				}
				if(pm_support > MAX_SUPPORT){ break;}
			}

			for(j=midCD;j<mp_reads->size;j++){
				if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->AB),vector_get(mp_reads,j),MAX_FRAG_SIZE)){
					mp_support++;				
				}
				if(SV_VECTOR_GET(dups,i)->AB.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
					break;
				}
				if(mp_support > MAX_SUPPORT){ break;}
			}
		}
		else if(SV_VECTOR_GET(dups,i)->orientation == DUP_BACK_COPY){

			midAB = scl_binary_search(pm_reads,&(SV_VECTOR_GET(dups,i)->AB));
			midCD = scl_binary_search(mp_reads,&(SV_VECTOR_GET(dups,i)->CD));
			for(j=midAB;j<pm_reads->size;j++){
				if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->AB),vector_get(pm_reads,j),MAX_FRAG_SIZE)){
					pm_support++;				
				}
				if(SV_VECTOR_GET(dups,i)->AB.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
					break;
				}
				if(pm_support > MAX_SUPPORT){ break;}	
			}

			for(j=midCD;j<mp_reads->size;j++){
				if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->CD),vector_get(mp_reads,j),MAX_FRAG_SIZE)){
					mp_support++;				
				}
				if(SV_VECTOR_GET(dups,i)->CD.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
					break;
				}

				if(mp_support > MAX_SUPPORT){ break;}	
			}
		}
		SV_VECTOR_GET(dups,i)->supports[0] = mp_support;
		SV_VECTOR_GET(dups,i)->supports[1] = pm_support;
	}

}

void update_inversion_supports_b(sv_t *inv, vector_t *pp_reads, vector_t *mm_reads){
	int j;
	int pp_support;
	int mm_support;
	int midAB;
	int midCD;


	pp_support = 0;
	mm_support = 0;
	midAB = scl_binary_search(pp_reads,&(inv->AB));
	midCD = scl_binary_search(mm_reads,&(inv->CD));


	for(j=midAB;j<pp_reads->size;j++){

		if(interval_pair_overlaps(&(inv->AB),vector_get(pp_reads,j),CLONE_MEAN/2)){
			pp_support++;				
		}
		if(inv->AB.end1 < IDIS_VECTOR_GET(pp_reads,j)->start1){
			break;
		}

		if(pp_support > MAX_SUPPORT){ break;}	
	}

	for(j=midCD;j<mm_reads->size;j++){
		if(interval_pair_overlaps(&(inv->CD),vector_get(mm_reads,j),CLONE_MEAN/2)){
			mm_support++;				
		}
		if(inv->CD.end1 < IDIS_VECTOR_GET(mm_reads,j)->start1){
			break;
		}

		if(mm_support > MAX_SUPPORT){ break;}	
	}

	inv->supports[0] = pp_support;
	inv->supports[1] = mm_support;


}


void update_inversion_supports(vector_t *inversions, vector_t *pp_reads, vector_t *mm_reads){
	int i,j;
	int pp_support;
	int mm_support;
	int midAB;
	int midCD;

	qsort(pp_reads->items,pp_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	qsort(mm_reads->items,mm_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	for(i=0;i<inversions->size;i++){
		pp_support = 0;
		mm_support = 0;
		midAB = scl_binary_search(pp_reads,&(SV_VECTOR_GET(inversions,i)->AB));
		midCD = scl_binary_search(mm_reads,&(SV_VECTOR_GET(inversions,i)->CD));


		for(j=midAB;j<pp_reads->size;j++){

			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->AB),vector_get(pp_reads,j),CLONE_MEAN)){
				pp_support++;				
			}
			if(SV_VECTOR_GET(inversions,i)->AB.end1 < IDIS_VECTOR_GET(pp_reads,j)->start1){
				break;
			}

			if(pp_support > MAX_SUPPORT){ break;}	
		}

		for(j=midCD;j<mm_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->CD),vector_get(mm_reads,j),CLONE_MEAN)){
				mm_support++;				
			}
			if(SV_VECTOR_GET(inversions,i)->CD.end1 < IDIS_VECTOR_GET(mm_reads,j)->start1){
				break;
			}

			if(mm_support > MAX_SUPPORT){ break;}	
		}

		SV_VECTOR_GET(inversions,i)->supports[0] = pp_support;
		SV_VECTOR_GET(inversions,i)->supports[1] = mm_support;

	}
}

void update_sv_supports_b(vector_t *svs, bam_vector_pack *reads){
	int i;

	qsort(reads->pp_discordants->items,reads->pp_discordants->size,sizeof(interval_discordant *),interval_pair_comp);
	qsort(reads->mm_discordants->items,reads->mm_discordants->size,sizeof(interval_discordant *),interval_pair_comp);
	qsort(reads->pm_discordants->items,reads->pm_discordants->size,sizeof(interval_discordant *),interval_pair_comp);
	qsort(reads->mp_discordants->items,reads->mp_discordants->size,sizeof(interval_discordant *),interval_pair_comp);

	for(i=0;i<svs->size;i++){
		sv_t *sv = vector_get(svs,i);
		switch(sv->type){
			case SV_INVERSION:
				update_inversion_supports_b(sv,reads->pp_discordants,reads->mm_discordants);
				break;
			case SV_DIRECT_DUPLICATION:
				update_duplication_supports_b(sv,reads->pm_discordants,reads->mp_discordants);
				break;
			case SV_INVERTED_DUPLICATION:
				update_duplication_supports_b(sv,reads->pp_discordants,reads->mm_discordants);
				break;
            case SV_TANDEM_DUPLICATION:
				update_tandem_duplication_supports_b(sv,reads->mp_discordants);
				break;
            case SV_TRANSLOCATION:
                update_duplication_supports_b(sv,reads->pm_discordants,reads->mp_discordants);
                break;
            case SV_INVERTED_TRANSLOCATION:
                update_duplication_supports_b(sv,reads->pp_discordants,reads->mm_discordants);
                break;
            case SV_DELETION:
				update_deletion_supports_b(sv,reads->pm_discordants);
				break;
			default:
				fprintf(stderr,"Unknown SV type ordinal %d\n",sv->type);
				VALOR_LOG("Unknown SV type ordinal %d\n",sv->type);
				exit(-1);
		}
	}

}

void update_sv_supports(vector_t *svs, bam_vector_pack *reads  ,sv_type type){

	switch(type){
		case SV_INVERSION:
			update_inversion_supports(svs,reads->pp_discordants,reads->mm_discordants);
			break;
		case SV_DIRECT_DUPLICATION:
			update_duplication_supports(svs,reads->pm_discordants,reads->mp_discordants);
			break;
		case SV_INVERTED_DUPLICATION:
			update_duplication_supports(svs,reads->pp_discordants,reads->mm_discordants);
			break;
		case SV_DELETION:
			update_deletion_supports(svs,reads->pm_discordants);
			break;
        case SV_TRANSLOCATION:
			update_duplication_supports(svs,reads->pm_discordants,reads->mp_discordants);
			break;
        case SV_INVERTED_TRANSLOCATION:
			update_duplication_supports(svs,reads->pp_discordants,reads->mm_discordants);
			break;
		default:
			fprintf(stderr,"Unknown SV type ordinal %d\n",type);
			VALOR_LOG("Unknown SV type ordinal %d\n",type);
			exit(-1);
	}

}
//TODO is this revers?
splitmolecule_t *inversion_reduce_breakpoints(sv_t *inv){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));

	sc->start1 = inv->AB.end1;
	sc->end1 = inv->CD.start1;
	sc->start2 = inv->AB.end2;
	sc->end2 = inv->CD.start2;
	return sc;
}
splitmolecule_t *duplication_reduce_breakpoints(sv_t *dup){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));
	if(dup->orientation == DUP_FORW_COPY){
		sc->start1 = dup->AB.start1;
		sc->end1 = dup->CD.end1;
		sc->start2 = dup->AB.end2;
		sc->end2 = dup->CD.start2;
	} else if(dup->orientation == DUP_BACK_COPY){

		sc->start1 = dup->AB.start2;
		sc->end1 = dup->CD.end2;
		sc->start2 = dup->AB.end1;;
		sc->end2 = dup->CD.start1;
	}
	return sc;
}

splitmolecule_t *inverted_duplication_reduce_breakpoints(sv_t *dup){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));
	if(dup->orientation == DUP_FORW_COPY){
		sc->start1 = dup->AB.start1;
		sc->end1 = dup->CD.end1;
		sc->start2 = dup->CD.end2;
		sc->end2 = dup->AB.start2;
	} else if(dup->orientation == DUP_BACK_COPY){
		sc->start2 = dup->AB.end1;
		sc->end2 = dup->CD.start1;
		sc->start1 = dup->CD.start2;
		sc->end1 = dup->AB.end2;
	}

	return sc;
}

splitmolecule_t *tandem_duplication_reduce_breakpoints(sv_t *sv){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));
	sc->start1 = sv->AB.start1;
	sc->end1 = sv->AB.end1;
	sc->start2 = sv->AB.start2;
    sc->end2 = sv->AB.end2;
	return sc;
}
splitmolecule_t *deletion_reduce_breakpoints(sv_t *sv){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));
	sc->start1 = sv->AB.start1;
	sc->start2 = sv->AB.start2;
	sc->end1 = sv->AB.end1;
	sc->end2 = sv->AB.end2;
	return sc;
}

splitmolecule_t *sv_reduce_breakpoints(sv_t *sv){
	switch(sv->type){
		case SV_INVERSION: return inversion_reduce_breakpoints(sv);

		case SV_DIRECT_DUPLICATION: return duplication_reduce_breakpoints(sv);
		case SV_INVERTED_DUPLICATION: return inverted_duplication_reduce_breakpoints(sv);
        case SV_TRANSLOCATION: return duplication_reduce_breakpoints(sv);
		case SV_INVERTED_TRANSLOCATION: return inverted_duplication_reduce_breakpoints(sv);
		case SV_DELETION: return deletion_reduce_breakpoints(sv);
		case SV_TANDEM_DUPLICATION: return tandem_duplication_reduce_breakpoints(sv);
		default: fprintf(stderr, "SV type of unknown ordinal %d!\n",sv->type); exit(-1);                                   ;
	}
}

int inversion_is_proper(sv_t *sv){
	sonic *snc = sonic_load(NULL);
	//bam_info *in_bams = get_bam_info(NULL);
	parameters *params = get_params();
	int chr = sv->chr;

    fprintf(logFile,"%s\t%d\t%d\t",snc->chromosome_names[chr],sv->AB.start1,sv->CD.end1);
    double ploidy = params->chr_copy_count[chr];
	if( params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],sv->AB.start1,sv->CD.end1)){
		fprintf(logFile,"sat 5'\n");
		return 0;
	}
	if( params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],sv->AB.start2,sv->CD.end2)){
	
		fprintf(logFile,"sat 3'\n");
		return 0;
	}
	if( sv->supports[0] < INVERSION_MIN_REQUIRED_SUPPORT / ploidy){
		fprintf(logFile,"sup 5'\n");
		return 0;
	}
	if( sv->supports[1] < INVERSION_MIN_REQUIRED_SUPPORT / ploidy){
		fprintf(logFile,"sup 3'\n");
		return 0;
	}

	
	if(sonic_is_gap(snc, snc->chromosome_names[chr], sv->AB.start1-CLONE_MEAN/2, sv->CD.end1+CLONE_MEAN/2) ||
					sonic_is_gap(snc, snc->chromosome_names[chr], sv->AB.start2-CLONE_MEAN/2, sv->CD.end2+CLONE_MEAN/2)){
			fprintf(logFile,"Gap\n");
			return 0;
		}

		if(sonic_is_gap(snc, snc->chromosome_names[chr], sv->AB.start1, sv->CD.end2)){
			int _start = sv->AB.start1;
			int _end = sv->CD.end2;
			sonic_interval *interval = sonic_intersect(snc,snc->chromosome_names[chr],_start,_end,SONIC_GAP);
			if(((double)interval->end-interval->start)/((double)_end-_start)>0.25){
			
				fprintf(logFile,"Gap%%\n");
				return 0;
		}

	}
	fprintf( logFile, "Call\n");
	return 1;
}
int invert_duplication_is_proper(sv_t *sv){
	parameters *params = get_params();
	



	sonic *snc = sonic_load(NULL);
	bam_info *in_bams = get_bam_info(NULL);

	int start = -1;
	int end = -1;
	int target_start = -1;
	int target_end = -1;
	int chr = sv->chr;	

	double ploidy = params->chr_copy_count[chr];
	
    if( sv->supports[0]<DUPLICATION_MIN_REQUIRED_SUPPORT/(ploidy)){
		return 0;
	}
	if( sv->supports[1]<DUPLICATION_MIN_REQUIRED_SUPPORT/(ploidy)){
		return 0;
	}
    if(sv->orientation ==DUP_FORW_COPY){
		target_start = sv->AB.start2;
		target_end = sv->CD.end2;
		start = sv->CD.start1;
		end = sv->AB.end1;
	}else if(sv->orientation == DUP_BACK_COPY){
		target_start = sv->AB.start1;
		target_end = sv->CD.end1;
		start = sv->CD.start2;
		end = sv->AB.end2;
	}	
	if(target_start > target_end){
//		int temp = target_start;
//		target_start = target_end;
//		target_end = temp;
	
        int avg = (target_start - target_end)/2 + target_end;

        target_start = avg - CLONE_MEAN /2;
        target_end = avg + CLONE_MEAN /2;
    }


	int is_ref_dup_source = sonic_is_segmental_duplication(snc,snc->chromosome_names[chr],start,end);
	int is_ref_dup_target = sonic_is_segmental_duplication(snc,snc->chromosome_names[chr],target_start,target_end);

	int is_ref_gap_source = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[chr],start,end);
	int is_ref_gap_target = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[chr],target_start,target_end);
	int is_ref_sat_source = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],start,end);
	int is_ref_sat_target = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],target_start,target_end);
	int does_cnv_support_dup;

//	does_cnv_support_dup= get_depth_region(in_bams->depths[chr],start,end) > (3.0/ploidy) * in_bams->depth_mean[chr] - (3.0/ploidy) * in_bams->depth_std[chr];

	does_cnv_support_dup= get_depth_region(in_bams->depths[chr],start,end) > in_bams->depth_mean[chr] + (3.0/ploidy) * in_bams->depth_std[chr];
	return !((is_ref_dup_source && is_ref_dup_target) || 
            !does_cnv_support_dup || 
            is_ref_gap_source || 
            is_ref_gap_target || 
            is_ref_sat_source || 
            is_ref_sat_target);
}


int direct_duplication_is_proper(sv_t *sv){
	parameters *params = get_params();

	sonic *snc = sonic_load(NULL);
	bam_info *in_bams = get_bam_info(NULL);

	int start = -1;
	int end = -1;
	int target_start = -1;
	int target_end = -1;
	int chr = sv->chr;

	double ploidy = params->chr_copy_count[chr];
	if( sv->supports[0]<DUPLICATION_MIN_REQUIRED_SUPPORT/(ploidy)){
		return 0;
	}
	if( sv->supports[1]<DUPLICATION_MIN_REQUIRED_SUPPORT/( ploidy)){
		return 0;
	}
    if(sv->orientation ==DUP_FORW_COPY){
        target_start = sv->AB.end2;
        target_end = sv->CD.start2;
        start = sv->AB.start1;
        end = sv->CD.end1;
    }else if(sv->orientation == DUP_BACK_COPY){
        target_start = sv->AB.end1;
        target_end = sv->CD.start1;
        start = sv->AB.start2;
        end = sv->CD.end2;
    }

	if(target_start > target_end){
//		int temp = target_start;
//		target_start = target_end;
//		target_end = temp;
        int avg = (target_start - target_end)/2 + target_end;

        target_start = avg - CLONE_MEAN /2;
        target_end = avg + CLONE_MEAN /2;
    }
	int is_ref_dup_source = sonic_is_segmental_duplication(snc,snc->chromosome_names[chr],start,end);
	int is_ref_dup_target = sonic_is_segmental_duplication(snc,snc->chromosome_names[chr],target_start,target_end);

	int is_ref_gap_source = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[chr],start,end);
	int is_ref_gap_target = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[chr],target_start,target_end);
	int is_ref_sat_source = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],start,end);
	int is_ref_sat_target = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[chr],target_start,target_end);

	int does_cnv_support_dup= get_depth_region(in_bams->depths[chr],start,end) >  in_bams->depth_mean[chr] + (3.0/ploidy) * in_bams->depth_std[chr];
//	int does_cnv_support_dup= get_depth_region(in_bams->depths[chr],start,end) > (3.0/ploidy) * in_bams->depth_mean[chr] - (3.0/ploidy) * in_bams->depth_std[chr];

    fprintf(logFile,"%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
            snc->chromosome_names[chr],start,end,
            snc->chromosome_names[chr],target_start,target_end,
            sv_type_name( sv->type));
    fprintf(logFile,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" ,get_depth_region(in_bams->depths[chr],start,end),does_cnv_support_dup,is_ref_dup_source, is_ref_dup_target, is_ref_gap_source, is_ref_gap_target, is_ref_sat_source, is_ref_sat_target);  
	return !((is_ref_dup_source && is_ref_dup_target) || 
            !does_cnv_support_dup || 
            is_ref_gap_source || 
            is_ref_gap_target || 
            is_ref_sat_source || 
            is_ref_sat_target);
}

int deletion_is_proper(sv_t *sv){

	sonic *snc = sonic_load(NULL);
	bam_info *in_bams = get_bam_info(NULL);
	parameters *params = get_params();
	int chr = sv->chr;
	double ploidy = params->chr_copy_count[chr];

	if( params->filter_gap && 
			sonic_is_gap(snc,snc->chromosome_names[chr], sv->AB.start1, sv->AB.end2)){
		return 0;
	}
	if( sv->supports[0] < DELETION_MIN_REQUIRED_SUPPORT/ploidy){
		return 0;
	}
	if( get_depth_region(in_bams->depths[chr],sv->AB.end1,sv->AB.start2) > (1 - 1.0 /ploidy) * in_bams->depth_mean[chr]  + (3.0 / ploidy) * in_bams->depth_std[chr]){
		return 0;
	}
	return 1;
}

int sv_is_proper(void *vsv){
	sv_t *sv = vsv;
	switch(sv->type){
	case SV_DELETION:
		return deletion_is_proper(sv);
	case SV_INVERSION:
		return inversion_is_proper(sv);
	case SV_DIRECT_DUPLICATION:
		return direct_duplication_is_proper(sv);
	case SV_INVERTED_DUPLICATION:
		return invert_duplication_is_proper(sv); 
	default:
		fprintf(stderr,"Sv type with ordinal %d is not implemented", (int)sv->type);
		exit(-1);
	}
}

