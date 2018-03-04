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
	default:
		return "unknown";
	}
}

void sv_fprint(FILE *stream, int chr, sv_t *t){
	 fprintf(stream,"chr:%d\t%d\t%d\t%d\t%d\t%lu\t-\t%d\t%d\t%d\t%d\t%lu\t++%d\t--%d\t%s\n",
                                chr,
                                t->AB.start1,
                                t->AB.end1,
                                t->AB.start2,
                                t->AB.end2,
                                t->AB.barcode,
                                t->CD.start1,
                                t->CD.end1,
                                t->CD.start2,
                                t->CD.end2,
                                t->CD.barcode,
				t->supports[0],
				t->supports[1],
				sv_type_name(t->type));
}

//TODO consider barcode
int sv_equals(const void *i1, const void *i2){
	const sv_t *ii1 = i1;
	const sv_t *ii2 = i2;
	return ii1->AB.start1 == ii2->AB.start1 &&
		ii1->AB.start2 == ii2->AB.start2 &&
		ii1->AB.end1 == ii2->AB.end1 &&
		ii1->AB.end2 == ii2->AB.end2 &&
		ii1->CD.start1 == ii2->CD.start1 &&
		ii1->CD.start2 == ii2->CD.start2 &&
		ii1->CD.end1 == ii2->CD.end1 &&
		ii1->CD.end2 == ii2->CD.end2 &&
		ii1->AB.barcode == ii2->AB.barcode &&
		ii1->CD.barcode == ii2->CD.barcode &&
		ii1->type == ii2->type;

//	return 	memcmp(&(ii1->AB),&(ii2->AB),sizeof(splitmolecule_t)) &&
//		memcmp(&(ii1->CD),&(ii2->CD),sizeof(splitmolecule_t));	
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

//TODO Sanity Check
int splitmolecule_indicates_duplication(splitmolecule_t *s1, splitmolecule_t *s2){
	
	if( !(s1->start1 < s2->start1 &&
		s1->end1 < s2->end1 &&
		s1->start2 < s2->start2 &&
		s1->end2 < s2->end2)){
		return 0;
	}
	if(     // S21-E11(
		i_distance(s2->start1, s1->start1, s2->end1, s1->end1) > DUP_OVERLAP &&
		i_distance(s2->start1, s1->start1, s2->end1, s1->end1) < DUP_GAP &&
		// E22-S12
		d_distance(s2->start1, s1->start2, s2->end1,s1->end2) < DUP_MAX_SIZE &&
		d_distance(s2->start1, s1->start2, s2->end1,s1->end2) > DUP_MIN_SIZE
	){
		return 1;
	}
	if(
		// S22-E12
		i_distance(s2->start2, s1->start2, s2->end1, s1->end2) > DUP_OVERLAP &&
		i_distance(s2->start2, s1->start2, s2->end1, s1->end2) < DUP_GAP &&
		// E21-S21
		d_distance(s1->start1, s2->start1, s2->end1,s1->end1) < DUP_MAX_SIZE &&
		d_distance(s1->start1, s2->start1, s2->end1,s1->end1) > DUP_MIN_SIZE
	){
		return 2;
	}
	return 0;
}

//TODO make sure different barcode
int splitmolecule_indicates_inversion(splitmolecule_t *s1, splitmolecule_t *s2){
	return s1->start1 < s2->start1 &&
		s1->end1 < s2->end1 &&
		s1->start2 < s2->start2 &&
		s1->end2 < s2->end2 &&
		i_distance(s2->start1,s1->start1,s2->end1,s1->end1) > INV_OVERLAP &&
		i_distance(s2->start2,s1->start2,s2->end2,s1->end2) > INV_OVERLAP &&
		i_distance(s2->start1,s1->start1,s2->end1,s1->end1) < INV_GAP &&
		i_distance(s2->start2,s1->start2,s2->end2,s1->end2) < INV_GAP;
}

int splitmolecule_indicates_sv(splitmolecule_t *s1, splitmolecule_t *s2, sv_type type){
	switch(type){
	case SV_INVERSION:
		return splitmolecule_indicates_inversion(s1,s2);
	case SV_DUPLICATION:
		return splitmolecule_indicates_duplication(s1,s2);
	default:
		fprintf(stderr,"Unknown SV type ordinal %d\n",type);
		VALOR_LOG("Unknown SV type ordinal %d\n",type);
		exit(-1);
	}
	return 0;
}

int TWO_ZEROS[2] = {0,0};
//double FOUR_ZEROS[4] = {0};
sv_t *sv_init(splitmolecule_t *sc1,splitmolecule_t *sc2,sv_type type){
	sv_t *new_i = getMem(sizeof(sv_t));
	memcpy(new_i->supports,TWO_ZEROS,sizeof(int)*2);
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
		});
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
		});
}
//Assumes i1 and i2 are same type SV
int sv_overlaps(sv_t *i1, sv_t *i2){
	switch(i1->type){
	case SV_INVERSION:
		return inversion_overlaps(i1,i2);
	case SV_DUPLICATION:
		return duplication_overlaps(i1,i2);
	default:
		fprintf(stderr,"Unknown SV type ordinal %d\n",i1->type);
		VALOR_LOG("Unknown SV type ordinal %d\n",i1->type);
		exit(-1);
	}
	return 0;
}
graph_t *make_sv_graph(vector_t *svs){
	graph_t *g = graph_init(svs->size * 2, sizeof(sv_t));
	g->hf = &sv_hf;

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
			ii = I10X_VECTOR_GET(regions,i);
			ij = I10X_VECTOR_GET(regions,j);
			if( interval_can_pair(ii,ij)){
				splitmolecule_t *new_molecule = splitmolecule_init(ii,ij);
				vector_soft_put(smolecules,new_molecule);
			}if( interval_can_pair(ij,ii)){
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
			if( splitmolecule_indicates_sv(AB,CD,type)){
				sv_t *tmp = sv_init(AB,CD,type);
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
		midAB = scl_binary_search(pm_reads,&(SV_VECTOR_GET(dups,i)->AB));
		midCD = scl_binary_search(mp_reads,&(SV_VECTOR_GET(dups,i)->CD));

		for(j=midAB;j<pm_reads->size;j++){
//			printf("midAB: %d, j %d",midAB,j);
			if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->AB),vector_get(pm_reads,j))){
				pm_support++;				
			}
//			printf("\n");
			if(SV_VECTOR_GET(dups,i)->AB.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
				break;
			}	
		}

		for(j=midCD;j<mp_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->CD),vector_get(mp_reads,j))){
				mp_support++;				
			}
			if(SV_VECTOR_GET(dups,i)->CD.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
				break;
			}	
		}
		midAB = scl_binary_search(pm_reads,&(SV_VECTOR_GET(dups,i)->CD));
		midCD = scl_binary_search(mp_reads,&(SV_VECTOR_GET(dups,i)->AB));

		for(j=midAB;j<pm_reads->size;j++){
//			printf("midAB: %d, j %d",midAB,j);
			if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->CD),vector_get(pm_reads,j))){
				pm_support++;				
			}
//			printf("\n");
			if(SV_VECTOR_GET(dups,i)->CD.end1 < IDIS_VECTOR_GET(pm_reads,j)->start1){
				break;
			}	
		}

		for(j=midCD;j<mp_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(dups,i)->AB),vector_get(mp_reads,j))){
				mp_support++;				
			}
			if(SV_VECTOR_GET(dups,i)->AB.end1 < IDIS_VECTOR_GET(mp_reads,j)->start1){
				break;
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
			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->AB),vector_get(pp_reads,j))){
				pp_support++;				
			}
//			printf("\n");
			if(SV_VECTOR_GET(inversions,i)->AB.end1 < IDIS_VECTOR_GET(pp_reads,j)->start1){
				break;
			}	
		}

		for(j=midCD;j<mm_reads->size;j++){
			if(interval_pair_overlaps(&(SV_VECTOR_GET(inversions,i)->CD),vector_get(mm_reads,j))){
				mm_support++;				
			}
			if(SV_VECTOR_GET(inversions,i)->CD.end1 < IDIS_VECTOR_GET(mm_reads,j)->start1){
				break;
			}	
		}

		SV_VECTOR_GET(inversions,i)->supports[0] = pp_support;
		SV_VECTOR_GET(inversions,i)->supports[1] = mm_support;

	}
}


void update_sv_supports(vector_t *svs, vector_t *support_vector1, vector_t *support_vector2,sv_type type){
	switch(type){
	case SV_INVERSION:
		update_inversion_supports(svs,support_vector1,support_vector2);
		break;
	case SV_DUPLICATION:
		update_duplication_supports(svs,support_vector1,support_vector2);
		break;
	default:
		fprintf(stderr,"Unknown SV type ordinal %d\n",type);
		VALOR_LOG("Unknown SV type ordinal %d\n",type);
		exit(-1);
	}

}
//TODO is this revers?
splitmolecule_t *sv_reduce_breakpoints(sv_t *inv){
	splitmolecule_t *sc = getMem(sizeof(splitmolecule_t));
	
	sc->start1 = inv->AB.end1;
	sc->end1 = inv->CD.start1;
	sc->start2 = inv->AB.end2;
	sc->end2 = inv->CD.start2;
	return sc;
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

