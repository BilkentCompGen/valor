#include "valorconfig.h"
#include "structural_variation.h"
#include <stdio.h>
#include "progress.h"


const char *sv_type_name(sv_type type){
	switch(type){
		case SV_INVERSION:
			return "inversion";
		case SV_DUPLICATION:
			return "duplication";
		case SV_INVERTED_DUPLICATION:
			return "inverted-duplication";
		default:
			return "unknown";
	}
}

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


void g_dfs_step(graph_t *g, vector_t *comp, sv_t *sv){
        sv->covered = 1;
        vector_put(comp,sv);
        vector_t *edges = graph_get_edges(g,sv);
        int i;
        for(i=0;i<edges->size;i++){
                sv_t **val = vector_get(edges,i);
                if(!(*val)->covered){
                        g_dfs_step(g,comp,*val);
                }
        }
}

vector_t *g_dfs_components(graph_t *g){
	adjlist_t *al = graph_to_al(g);
        vector_t *comps = vector_init(sizeof(vector_t),16);
	comps->rmv  = &vector_free;
        int i;

        for(i=0;i<al->size;i++){
		sv_t *sv = al_get_value(al,i);
                if( !sv->covered){
                        vector_t *comp = vector_init(sizeof(sv_t),40);
                        comp->rmv = &sv_destroy;
			g_dfs_step(g,comp,sv);
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

//Distance for inversions
int i_distance(int start1, int start2, int end1, int end2){
	if(start1 < start2){
		return i_distance(start2,start1,end2,end1);
	}
	return start1-end2;
}

//Distance for duplications
int d_distance(int start1, int start2, int end1, int end2){
	if(start1 < start2){
		return d_distance(start2, start1, end2, end1);
	}
	return end1-start2;
}

int interval_outer_distance(interval_10X a, interval_10X b){
	if( a.end > b.end){
		return a.end - b.start;
	}
	return b.end - a.start;
}
int interval_inner_distance(interval_10X a, interval_10X b){
	if( a.start < b.start){
		return b.start - a.end;
	}
	return a.start - b.end;
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
		case SV_DUPLICATION:
			return splitmolecule_indicates_duplication(*s1,*s2);
		case SV_INVERTED_DUPLICATION:
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

	//	new_i->depths[0] = 0;
	//	new_i->depths[1] = 1;
	//	new_i->depths[2] = 2;
	//	new_i->depths[3] = 3;
	if(sc2!=NULL){//In case I use this function for sv's with 1 split molecule
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


splitmolecule_t *splitmolecule_init(interval_10X *i1, interval_10X *i2){
	if( i1->barcode != i2->barcode){
		fprintf(stderr,"Can't create split molecule of different barcodes\n");
		return NULL;
	}
	splitmolecule_t *new_molecule = malloc(sizeof(splitmolecule_t));
	new_molecule->start1 = i1->start;
	new_molecule->start2 = i2->start;
	new_molecule->end1 = i1->end;
	new_molecule->end2 = i2->end;
	new_molecule->barcode = i1->barcode;
	return new_molecule;
}

void splitmolecule_destroy(splitmolecule_t *molecule){
	free(molecule);
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
//Assumes i1 and i2 are same type SV
int sv_overlaps(sv_t *i1, sv_t *i2){
	switch(i1->type){
		case SV_INVERSION:
			return inversion_overlaps(i1,i2);
		case SV_DUPLICATION:
			return duplication_overlaps(i1,i2);
		case SV_INVERTED_DUPLICATION:
			return duplication_overlaps(i1,i2);
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
//                        if(!graph_have_node(g,*ptr)){
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
	printf("Adding Nodes\n");
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


#define SCL_INIT_LIMIT 16


vector_t *discover_split_molecules(vector_t *regions){
	int i,j;
	vector_t *smolecules = vector_init(sizeof(splitmolecule_t),SCL_INIT_LIMIT);
	interval_10X *ii;
	interval_10X *ij;
	printf("Discoverng Split Molecules...\n\n");
	for(i=0;i<regions->size;i++){
		for(j=i+1;j<regions->size;j++){
			ii = vector_get(regions,i);
			ij = vector_get(regions,j);
			if( interval_can_pair(ii,ij)){
				splitmolecule_t *new_molecule = splitmolecule_init(ii,ij);
				vector_soft_put(smolecules,new_molecule);
			}else if( interval_can_pair(ij,ii)){
				splitmolecule_t *new_molecule = splitmolecule_init(ij,ii);
				vector_soft_put(smolecules,new_molecule);
			}
			//	VALOR_LOG("%d-%d\t%d-%d\n",i,j,ii->barcode,ij->barcode);
			if( ij->barcode != ii->barcode){
				break;
			}	
		}
		if((i&511) == 0){
			update_progress(i,regions->size);
		}
	}

	return smolecules;
}

vector_t *find_svs(vector_t *split_molecules, sv_type type){
	int i,j;
	vector_t *svs = vector_init(sizeof(sv_t),SCL_INIT_LIMIT);

	splitmolecule_t *AB;
	splitmolecule_t *CD;

	printf("Discovering SV's of type %s\n\n", sv_type_name(type));

	for(i=0;i<split_molecules->size;i++){
		AB = vector_get(split_molecules,i);
		for(j=0;j<split_molecules->size;j++){
			if(i==j){continue;}
			CD = vector_get(split_molecules,j);
			char orient = splitmolecule_indicates_sv(AB,CD,type);
			if(orient){
				sv_t *tmp = sv_init(AB,CD,type);
				tmp->orientation = orient;
				vector_soft_put(svs,tmp);
			}
		}
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


size_t scl_binary_search(vector_t *intervals, splitmolecule_t *key){
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

void update_duplication_supports(vector_t *dups, vector_t *pm_reads, vector_t *mp_reads){
	int i,j;
	int pm_support;
	int mp_support;
	int midAB;
	int midCD;
	//	for(i=0;i<pp_reads->size;i++){
	//		interval_discordant *d = vector_get(pp_reads,i);
	//		printf("%d\t%d\t%d\t%d\n",d->start1,d->end1,d->start2,d->end2);
	//	}

	qsort(pm_reads->items,pm_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	//	for(i=0;i<pp_reads->size;i++){
	//		interval_discordant *d = vector_get(pp_reads,i);
	//		printf("%d\t%d\t%d\t%d\n",d->start1,d->end1,d->start2,d->end2);
	//	}

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

void update_inversion_supports(vector_t *inversions, vector_t *pp_reads, vector_t *mm_reads){
	int i,j;
	int pp_support;
	int mm_support;
	int midAB;
	int midCD;
	//	for(i=0;i<pp_reads->size;i++){
	//		interval_discordant *d = vector_get(pp_reads,i);
	//		printf("%d\t%d\t%d\t%d\n",d->start1,d->end1,d->start2,d->end2);
	//	}

	qsort(pp_reads->items,pp_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	//	for(i=0;i<pp_reads->size;i++){
	//		interval_discordant *d = vector_get(pp_reads,i);
	//		printf("%d\t%d\t%d\t%d\n",d->start1,d->end1,d->start2,d->end2);
	//	}

	qsort(mm_reads->items,mm_reads->size,sizeof(interval_discordant *),interval_pair_comp);

	for(i=0;i<inversions->size;i++){
		pp_support = 0;
		mm_support = 0;
		midAB = scl_binary_search(pp_reads,&(SV_VECTOR_GET(inversions,i)->AB));
		midCD = scl_binary_search(mm_reads,&(SV_VECTOR_GET(inversions,i)->CD));


		for(j=midAB;j<pp_reads->size;j++){
			//			printf("midAB: %d, j %d",midAB,j);
			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->AB),vector_get(pp_reads,j),MAX_FRAG_SIZE)){
				pp_support++;				
			}
			//			printf("\n");
			if(SV_VECTOR_GET(inversions,i)->AB.end1 < IDIS_VECTOR_GET(pp_reads,j)->start1){
				break;
			}
		
			if(pp_support > MAX_SUPPORT){ break;}	
		}

		for(j=midCD;j<mm_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->CD),vector_get(mm_reads,j),MAX_FRAG_SIZE)){
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


void update_sv_supports(vector_t *svs, bam_vector_pack *reads  ,sv_type type){
	switch(type){
		case SV_INVERSION:
			update_inversion_supports(svs,reads->pp_discordants,reads->mm_discordants);
			break;
		case SV_DUPLICATION:
			update_duplication_supports(svs,reads->pm_discordants,reads->mp_discordants);
			break;
		case SV_INVERTED_DUPLICATION:
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
splitmolecule_t *sv_reduce_breakpoints(sv_t *sv){
	switch(sv->type){
		case SV_INVERSION: return inversion_reduce_breakpoints(sv);
		case SV_DUPLICATION: return duplication_reduce_breakpoints(sv);
		case SV_INVERTED_DUPLICATION: return inverted_duplication_reduce_breakpoints(sv);
		default: return NULL;
	}
}
splitmolecule_t *splitmolecule_copy(splitmolecule_t *scl){
	splitmolecule_t *ncl = getMem(sizeof(splitmolecule_t));

	ncl->start1 = scl->start1;
	ncl->start2 = scl->start2;
	ncl->end1 = scl->end1;
	ncl->end2 = scl->end2;
	ncl->barcode = scl->barcode;
	return ncl;
}

