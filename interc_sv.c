#include "interc_sv.h"
#include "graph.h"
#include "clique_inter.h"
#define SCL_INIT_LIMIT 400
#include "cnv.h"

int isv_chr_cmp(const ic_sv_t *s1, const ic_sv_t *s2){
    if( s1->chr_source == s2->chr_source){
        return s1->chr_target - s2->chr_target;
    }
    return s1->chr_source - s2->chr_source;
}

int isv_pos_cmp(const ic_sv_t *s1, const ic_sv_t *s2){
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
    if(s1->EF.start1 != s2->EF.start1){
        return s1->EF.start1 - s2->EF.start1;
    }
    if(s1->EF.end1 != s2->EF.end1){
        return s1->EF.end1 - s2->EF.end1;
    }

    if(s1->EF.start2 != s2->EF.start2){
        return s1->EF.start2 - s2->EF.start2;
    }
    if(s1->EF.end2 != s2->EF.end2){
        return s1->EF.end2 - s2->EF.end2;
    }

	return 0;
}
int _ic_sv_cmp(const void *v1, const void *v2, size_t size){
    const ic_sv_t *s1 = (void *)v1;
    const ic_sv_t *s2 = (void *)v2;
    int chr_dif = isv_chr_cmp(s1,s2);
    if(chr_dif != 0 ){
        return chr_dif;
    }
    int pos_dif = isv_pos_cmp(s1,s2);
    if(pos_dif != 0){
        return pos_dif;
    }

    return 0;
}


int interc_sv_comp(const void *v1, const void *v2){
    const ic_sv_t *s1 = *(void **)v1;
    const ic_sv_t *s2 = *(void **)v2;
    int chr_dif = isv_chr_cmp(s1,s2);
    if(chr_dif != 0 ){
        return chr_dif;
    }
    int pos_dif = isv_pos_cmp(s1,s2);
    if(pos_dif != 0){
        return pos_dif;
    }

    return 0;
}

int interc_sv_compd(const void *v1, const void *v2, size_t size){
    return interc_sv_comp(v1,v2);
}
size_t split_molecule_binary_search(vector_t *splits, interval_10X key){
	if(splits->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = splits->size - 1;

	while( first < last){
		mid = (first + last)/2;
		if(((interval_pair*)vector_get(splits,mid))->start1< key.start){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
	}
//-30000
	while( mid >= 0 && key.start -30000 < ((interval_pair*) vector_get(splits,mid))->start1){
		mid--;
	}
	return mid +1;
}

// So that program can safely assume no NULL vectors.
void filter_dangling_reads(vector_t *reads){
	int i;
	reads->REMOVE_POLICY = REMP_LAZY;
	bit_set_t *bs = get_bam_info(NULL)->chro_bs;
	for(i=0;i<reads->size;i++){
		barcoded_read_pair *p = vector_get(reads,i);

		if(bit_set_get_bit(bs,p->r_chr) == 0){
			vector_remove(reads,i);
		}
	}
	vector_defragment(reads);
	reads->REMOVE_POLICY = REMP_SORTED;
}
size_t interval_pair_binary_search(vector_t *mols, interval_10X interval){
	if(mols->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = mols->size - 1;

	while( first < last){
//		mid = (first + last)/2;
        mid = (last - first)/2 + first;
		if(((interval_pair*)vector_get(mols,mid))->start1 < interval.start){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
	}
    int cur_start = ((interval_pair *) vector_get(mols,mid))->start1;
	mid--;
	while( mid >= 0 && cur_start > interval.start){
		mid--;
	}
	return mid+1;
}

size_t molecule_barcode_binary_search(vector_t *mols, unsigned long key){
	if(mols->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = mols->size - 1;

	while( first < last){
//		mid = (first + last)/2;
        mid = (last - first)/2 + first;
		if(((interval_10X*)vector_get(mols,mid))->barcode < key){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
	}
	unsigned long cur_barcode = ((interval_10X *) vector_get(mols,mid))->barcode;
	mid--;
	while( mid >= 0 && cur_barcode == ((interval_10X *) vector_get(mols,mid))->barcode){
		mid--;
	}
	return mid+1;
}
void inter_sv_call_bed_print(FILE *stream, inter_sv_call_t *sv){
	sonic *snc = sonic_load(NULL);

    fprintf(stream,"%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n",
            snc->chromosome_names[sv->break_points.chr1],
            sv->break_points.start1,
            sv->break_points.end1,
            snc->chromosome_names[sv->break_points.chr2],
            sv->break_points.start2,
            sv->break_points.end2,
            sv_type_name(sv->type),
            sv->cluster_size,
            sv->supports[0],
            sv->supports[1],
            sv->supports[2]
           );
}

void ic_sv_bed_print(FILE *stream, ic_sv_t *sv){
	sonic *snc = sonic_load(NULL);
    if( sv->chr_target == sv->AB.chr1){
        fprintf(stream,"%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
                snc->chromosome_names[sv->chr_source],
                sv->EF.end1,
                sv->EF.start2,
                snc->chromosome_names[sv->chr_target],
                sv->AB.end1,
                sv->CD.start1,
                sv_type_name(sv->type)
               );
    }else{

        fprintf(stream,"%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
                snc->chromosome_names[sv->chr_source],
                sv->EF.end1,
                sv->EF.start2,
                snc->chromosome_names[sv->chr_target],
                sv->AB.end2,
                sv->CD.start1,
                sv_type_name(sv->type)
               );
    }
}
int inter_split_overlaps(inter_split_molecule_t s1, inter_split_molecule_t s2, int relax){

    if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}

    return interval_overlaps(&(interval_10X){s1.start1,s1.end1,0},&(interval_10X){s2.start1,s2.end1,0},relax) &&
            interval_overlaps(&(interval_10X){s1.start2,s1.end2,0},&(interval_10X){s2.start2,s2.end2,0},relax);

}


int direct_translocation_overlaps(ic_sv_t *i1, ic_sv_t *i2){
    return  i1->chr_source == i2->chr_source && i1->chr_target == i2->chr_target &&
        interval_pair_overlaps(&i1->EF,&i2->EF,2*CLONE_MEAN) && 
        inter_split_overlaps(i1->AB,i2->AB,2*CLONE_MEAN) &&
        inter_split_overlaps(i1->CD,i2->CD,2*CLONE_MEAN);
    return 0;
}
 
int invert_translocation_overlaps(ic_sv_t *i1, ic_sv_t *i2){
    return 0;
}
 

int inter_sv_overlaps(ic_sv_t *i1, ic_sv_t *i2){

	if(i1->type != i2->type) { return 0;}
    switch(i1->type){
        case SV_TRANSLOCATION:
            return direct_translocation_overlaps(i1,i2);
        case SV_INVERTED_TRANSLOCATION:
            return invert_translocation_overlaps(i1,i2);
        default:
            fprintf(stderr,"Invalid Inter SV ordinal: %d\n",i1->type);
            exit(-1);
    }
    return 0;
}
inter_split_molecule_t *inter_split_init(barcoded_read_pair *pair, interval_10X *a, interval_10X *b){
	inter_split_molecule_t *isplit = malloc(sizeof(inter_split_molecule_t));

	isplit->start1 = a->start;
	isplit->start2 = b->start;
	isplit->end1 = a->end;
	isplit->end2 = b->end;
	isplit->chr1 = pair->l_chr;
	isplit->chr2 = pair->r_chr;
    isplit->barcode = pair->barcode;
	return isplit;
}

int inter_split_indicates_invert_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2){
	if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}
    if( s1.barcode == s2.barcode) {return 0;}

	interval_10X A = (interval_10X){s1.start1,s1.end1,s1.barcode};  
	interval_10X B = (interval_10X){s1.start2,s1.end2,s1.barcode};  
	interval_10X C = (interval_10X){s2.start1,s2.end1,s2.barcode};  
	interval_10X D = (interval_10X){s2.start2,s2.end2,s2.barcode};  
	if(interval_inner_distance(A,C) < DUP_GAP &&
			interval_inner_distance(A,C) > DUP_OVERLAP &&
            interval_outer_distance(B,D) < DUP_MAX_DIST){
		return INTERC_BACK_COPY;
	}
	else if(interval_inner_distance(B,D) < DUP_GAP &&
			interval_inner_distance(B,D) > DUP_OVERLAP &&
              interval_outer_distance(A,C) < DUP_MAX_DIST){
		return INTERC_FORW_COPY;
	}
	return 0; 
}
int inter_split_indicates_direct_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2){
	if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}
    if( s1.barcode == s2.barcode) {return 0;}

	interval_10X A = (interval_10X){s1.start1,s1.end1,s1.barcode};  
	interval_10X B = (interval_10X){s1.start2,s1.end2,s1.barcode};  
	interval_10X C = (interval_10X){s2.start1,s2.end1,s2.barcode};  
	interval_10X D = (interval_10X){s2.start2,s2.end2,s2.barcode};  
	if(interval_inner_distance(A,C) < DUP_GAP &&
			interval_inner_distance(A,C) > DUP_OVERLAP &&
            interval_outer_distance(B,D) < DUP_MAX_DIST){
		return INTERC_BACK_COPY;
	}
	else if(interval_inner_distance(B,D) < DUP_GAP &&
			interval_inner_distance(B,D) > DUP_OVERLAP &&
              interval_outer_distance(A,C) < DUP_MAX_DIST){
		return INTERC_FORW_COPY;
	}
	return 0; 
}


int inter_split_indicates_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2, sv_type type){
    switch(type){
    case SV_TRANSLOCATION:
        return inter_split_indicates_direct_translocation(s1,s2);
    case SV_INVERTED_TRANSLOCATION:
        return inter_split_indicates_invert_translocation(s1,s2);

    default:
        fprintf(stderr,"Inter SV with  unknown ordinal: %d!\n",type);
        exit(-1);
    }
    return 0;
}

ic_sv_t *inter_sv_init(inter_split_molecule_t *a, inter_split_molecule_t *b, interval_pair *tra_del, sv_type type, int orient){
	if(orient == 0){ return NULL;}
    ic_sv_t *new_i = malloc(sizeof(ic_sv_t));
	memset(new_i,0,sizeof(ic_sv_t));
	new_i->supports[0] = 1;//TODO fix this
	new_i->supports[1] = 1;//TODO fix this
	new_i->supports[2] = 1;//TODO fix this

	new_i->covered = 0;
	new_i->tabu = 0;
	new_i->dv = 0;
	new_i->inactive = 0;
	new_i->AB=*a;
	new_i->CD=*b;
	new_i->EF=*tra_del;
	new_i->type=type;


    if( orient == INTERC_FORW_COPY){
        new_i->chr_source = MIN(a->chr1,a->chr2);
        new_i->chr_target = MAX(a->chr1,a->chr2);
    }else if(orient == INTERC_BACK_COPY){
        new_i->chr_target = MIN(a->chr1,a->chr2);
        new_i->chr_source = MAX(a->chr1,a->chr2);

    }


	return new_i;
}

interval_pair *find_matching_split_molecule(vector_t **splits,inter_split_molecule_t *a, inter_split_molecule_t *b, int orient){     
	if(orient == 0){
       
    //    VALOR_LOG("BAD Orient: %d\n",orient);
        return NULL;}
    int src_chr; 
    int tar_chr;

    interval_10X deletion_interval = {0,0,0};
    if( orient == INTERC_FORW_COPY){
        src_chr = MIN(a->chr1,a->chr2);
        tar_chr = MAX(a->chr1,a->chr2);
        int start = MIN(a->end1,b->end1);
        int end = MAX(a->start1,b->start1);
        deletion_interval = (interval_10X){.start = start,.end = end,.barcode=0};
    }else if(orient == INTERC_BACK_COPY){
        tar_chr = MIN(a->chr1,a->chr2);
        src_chr = MAX(a->chr1,a->chr2);

        int start = MIN(a->end2,b->end2);
        int end = MAX(a->start2,b->start2);
        deletion_interval = (interval_10X){.start = start,.end = end,.barcode=0};
    }else{
        
        fprintf(stderr, "Translocation orientation cannot be %d, terminating!!\n",orient);
        exit(-1);
        return NULL;
    }
    //interval_pair *aaa = splitmolecule_init(&(interval_10X){deletion_interval.start-20000,deletion_interval.start,0},
     //                                       &(interval_10X){deletion_interval.end,deletion_interval.end+20000,0});
    //return aaa;
    size_t pos = split_molecule_binary_search(splits[src_chr],deletion_interval);
    if( pos == -1){
  //      VALOR_LOG("NO %d %d\n",deletion_interval.start,deletion_interval.end);
        return NULL;
    }

    interval_pair *cand = vector_get(splits[src_chr],pos);
    while( pos < splits[src_chr]->size && cand->start1 < deletion_interval.end + 50000){


        if( cand->start1 + cand->end1 < 2 * deletion_interval.start  && 
             cand->start1 + cand->end1 > 2 * deletion_interval.start - CLONE_MEAN  && 
                cand->start2 + cand->end2 < 2 * deletion_interval.end + CLONE_MEAN &&
                cand->start2 + cand->end2 > 2 * deletion_interval.end){
    
            return vector_get(splits[src_chr], pos);
        }
        pos++;
        cand = vector_get(splits[src_chr],pos);
    }

//    VALOR_LOG("CM\n%d\t%d\t%d\n",src_chr,deletion_interval.start,deletion_interval.end);
    return NULL;
}

vector_t *find_interc_translocations(vector_t *sp1, vector_t *sp2, vector_t **molecules,sv_type type){
	vector_t *tlocs = vector_init(sizeof(ic_sv_t),512);
	int i,j;

	for(i=0;i<sp1->size;i++){
		inter_split_molecule_t *a = vector_get(sp1,i);
		for(j=0;j<sp2->size;j++){   
			inter_split_molecule_t *b = vector_get(sp2,j);
			int orient = inter_split_indicates_translocation(*a,*b,type);

            interval_pair *del_tra = find_matching_split_molecule(molecules,a,b,orient);
			
//            VALOR_LOG("%d\t%d\n",orient,del_tra!=NULL);
            if(orient && del_tra != NULL){
				vector_soft_put(tlocs,inter_sv_init(a,b,del_tra,type,orient));
			}
		}
	}

	return tlocs;
}

int is_inter_chr_split(barcoded_read_pair *pair, interval_10X *a, interval_10X *b){
	if( (pair->barcode !=a->barcode) || (pair->barcode!=b->barcode)){ return 0;}

	return (in_range(pair->left,a->start,2*MOLECULE_EXT) && in_range(pair->right,b->start,2*MOLECULE_EXT));        
}

size_t ip_binary_search(vector_t *intervals, interval_pair *key){
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
    while(mid>0 && key->end1 > IDIS_VECTOR_GET(intervals,mid)->start1){
        mid--;
    }
	return mid;
}

int split_get_pm_support(interval_pair *split, vector_t *discordants){
    int support = 0;
    int mid = ip_binary_search (discordants,split);
    int j;
    for(j=mid;j < discordants->size;j++){
		if(interval_pair_overlaps(split,vector_get(discordants,j),CLONE_MEAN/2)){
			support++;				
		}
		if(split->end1 < IDIS_VECTOR_GET(discordants,j)->start1 + CLONE_MEAN/2){
			break;
		}
		if(support > MAX_SUPPORT){ break;}	
	}
    return support;
}

void filter_unsupported_pm_splits(vector_t *splits, vector_t *discordants){
    int i;


    splits->REMOVE_POLICY = REMP_LAZY;
    for(i=0;i<splits->size;i++){
        interval_pair *split = vector_get(splits,i);
        int support = split_get_pm_support(split,discordants);
        if( support < 1){
            vector_remove(splits,i);
        }
    } 
    splits->REMOVE_POLICY = REMP_SORTED;
    vector_defragment(splits);
    vector_zip(splits);    
}


// returns a vector of inter_split_molecules
vector_t *find_inter_split_molecules(vector_t *reads, vector_t *mol_a, vector_t *mol_b){
	vector_t *isms = vector_init(sizeof(inter_split_molecule_t),512);
	int i;

	for(i=0;i<reads->size;i++){
		barcoded_read_pair *pair =  vector_get(reads,i);
		if( pair->barcode == -1){ continue;}
        
        
/*
        size_t k = 0;

        while (k < mol_a->size && ((interval_10X *)vector_get(mol_a,k))->barcode != pair->barcode){
            k++;
        }
        if ( k >= mol_a->size){
            continue;
        }
        size_t t = 0;
        while (t < mol_b->size && ((interval_10X *)vector_get(mol_b,t))->barcode != pair->barcode){
            t++;
        }
        if ( t >= mol_b->size){
            continue;
        }

*/
		size_t k = molecule_barcode_binary_search(mol_a,pair->barcode);
        size_t t = molecule_barcode_binary_search(mol_b,pair->barcode);
        size_t kk = k;
		size_t tt = t;

		while(I10X_VECTOR_GET(mol_a,kk)->barcode <= pair->barcode){	
			while(I10X_VECTOR_GET(mol_b,tt)->barcode <= pair->barcode){

                if(I10X_VECTOR_GET(mol_a,kk)->barcode != pair->barcode ||
                    I10X_VECTOR_GET(mol_b,tt)->barcode != pair->barcode){
   
				    tt++;
				    if(tt >=mol_b->size){ break;}	
                    continue;
                } 
                        
                if(is_inter_chr_split(pair,vector_get(mol_a,kk),vector_get(mol_b,tt))){
                    vector_soft_put(isms,inter_split_init(pair,vector_get(mol_a,kk),vector_get(mol_b,tt)));
			    }else  if(is_inter_chr_split(pair,vector_get(mol_b,tt),vector_get(mol_a,kk))){
                    vector_soft_put(isms,inter_split_init(pair,vector_get(mol_b,tt),vector_get(mol_a,kk)));
				}	

                
				tt++;
				if(tt >=mol_b->size){ break;}	
			}
			tt = t;
			kk++;
			if(kk >= mol_a->size){ break;}
		}
	}
	return isms; 
}



void ic_sv_g_dfs_step(graph_t *g, vector_t *comp, ic_sv_t *sv){
	sv->covered = 1;
	vector_put(comp,sv);
	vector_t *edges = graph_get_edges(g,sv);
	int i;
	for(i=0;i<edges->size;i++){
		ic_sv_t **val = vector_get(edges,i);
		if(!(*val)->covered){
			ic_sv_g_dfs_step(g,comp,*val);
		}
	}
}
#define BIG_PRIME 1300501
size_t ic_sv_hf(hashtable_t *table, const void *vsv){
	const ic_sv_t *sv = vsv;
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


	//return SuperFastHash((vsv),2*sizeof(interval_pair)) % table->size;
}

vector_t *ic_sv_g_dfs_components(graph_t *g){
	adjlist_t *al = graph_to_al(g);
	vector_t *comps = vector_init(sizeof(vector_t),16);
	comps->rmv  = &vector_free;
	int i;

	for(i=0;i<al->size;i++){
		ic_sv_t *sv = al_get_value(al,i);
		if( !sv->covered){
			vector_t *comp = vector_init(sizeof(ic_sv_t),40);

			ic_sv_g_dfs_step(g,comp,sv);
			vector_soft_put(comps,comp);
		}

	}
	vector_free(al);
	return comps;
}
inter_interval_pair ic_sv_reduce_breakpoints(ic_sv_t *sv){
    if(sv->chr_target == sv->AB.chr1){
        return (inter_interval_pair){.chr1= sv->chr_source,
                                        .chr2=sv->chr_target,
                                        .start1=sv->EF.end1,
                                        .end1=sv->EF.start2,
                                        .start2=sv->AB.end1,
                                        .end2= sv->CD.start1,
                                        .barcode= 0};
    }else{
        return (inter_interval_pair){.chr1= sv->chr_source,
                                        .chr2=sv->chr_target,
                                        .start1=sv->EF.end1,
                                        .end1=sv->EF.start2,
                                        .start2=sv->AB.end2,
                                        .end2= sv->CD.start1,
                                        .barcode= 0};
    }
}

int ic_sv_is_proper(void *vcall){
    inter_sv_call_t *call =vcall;
    bam_info *in_bams = get_bam_info(NULL);
    parameters *params = get_params();
	sonic *snc = sonic_load(NULL);

    int start = call->break_points.start1;
    int target_start = call->break_points.start2;
    int end = call->break_points.end1;
    int target_end = call->break_points.end2;

    int src_chr = call->break_points.chr1;
    int tgt_chr = call->break_points.chr2;

    int is_ref_dup_source = sonic_is_segmental_duplication(snc,snc->chromosome_names[src_chr],start,end);
	int is_ref_dup_target = sonic_is_segmental_duplication(snc,snc->chromosome_names[tgt_chr],target_start,target_end);

	int is_ref_gap_source = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[src_chr],start,end);
	int is_ref_gap_target = params->filter_gap && sonic_is_gap(snc,snc->chromosome_names[tgt_chr],target_start,target_end);
	int is_ref_sat_source = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[src_chr],start,end);
	int is_ref_sat_target = params->filter_satellite && sonic_is_satellite(snc,snc->chromosome_names[tgt_chr],target_start,target_end);
    double depth = get_depth_region(in_bams->depths[src_chr],start,end);

    int does_cnv_support_tra;  
    if(is_ref_dup_source){
        does_cnv_support_tra=  get_depth_region(in_bams->depths[src_chr],start,end) > in_bams->depth_mean[src_chr] - 1.5 * in_bams->depth_std[src_chr];
   
    }else{

        does_cnv_support_tra=  depth < in_bams->depth_mean[src_chr] + 1.5 * in_bams->depth_std[src_chr] &&
            get_depth_region(in_bams->depths[src_chr],start,end) > in_bams->depth_mean[src_chr] - 1.5 * in_bams->depth_std[src_chr];
   
    }

    inter_sv_call_bed_print(logFile,call);  
    fprintf(logFile,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" ,depth,does_cnv_support_tra,is_ref_dup_source, is_ref_dup_target, is_ref_gap_source, is_ref_gap_target, is_ref_sat_source, is_ref_sat_target);  
    return !(is_ref_dup_source && is_ref_dup_target) && !(is_ref_gap_source || is_ref_gap_target) && !(is_ref_sat_source && is_ref_sat_target) && does_cnv_support_tra;
}
inter_sv_call_t *ic_sv_cluster_resolve(vector_t *cluster){
    inter_interval_pair bp;
    int supports[3] = {0};
    int i;
    int j;
    sv_type type = -1;
    int cnt = 0;
    for(i = 0;i< cluster->size;i++){
        ic_sv_t *sv =*(ic_sv_t **) vector_get(cluster,i);

        inter_interval_pair sv_bp =  ic_sv_reduce_breakpoints(sv);

        if( cnt != 0 && !inter_split_overlaps(bp,sv_bp,CLONE_MEAN)){
        //    continue;
        }
        if(cnt == 0){
            bp.chr1 = sv->chr_source;
            bp.chr2 = sv->chr_target;
            type = sv->type;
        }
         
        cnt++;

        bp.start1 = (double)(bp.start1)*((double)(cnt-1)/cnt) + (double)sv_bp.start1/cnt;
        bp.start2 = (double)(bp.start2)*((double)(cnt-1)/cnt) + (double)sv_bp.start2/cnt;
        bp.end1 = (double)(bp.end1)*((double)(cnt-1)/cnt) + (double)sv_bp.end1/cnt;
        bp.end2 = (double)(bp.end2)*((double)(cnt-1)/cnt) + (double)sv_bp.end2/cnt;
        for(j=0;j<3;j++){
            supports[j]+= sv->supports[j];
        }   
    }
    inter_sv_call_t *call = malloc(sizeof(inter_sv_call_t));
    call->break_points = bp;
    
    for(j=0;j<3;j++){
        call->supports[j] = supports[j];
    }
    call->cluster_size = cnt;
    call->type = type;
    return call;
}

vector_t *cluster_interchromosomal_events(vector_t *predictions){

    vector_t *calls = vector_init(sizeof(inter_sv_call_t),predictions->size);
    graph_t *sv_graph = graph_init(predictions->size *2, sizeof(ic_sv_t));
    sv_graph->hf = &ic_sv_hf;
    sv_graph->key_cmp = &_ic_sv_cmp; 

	parameters *params = get_params();
    int i,j;
    for(i=0; i<predictions->size; i++){
        graph_put_node(sv_graph,vector_get(predictions,i));
    }

    for(i=0; i<predictions->size; i++){
        ic_sv_t *a = vector_get(predictions,i);
        for(j=i+1; j<predictions->size; j++){
            ic_sv_t *b = vector_get(predictions,j);
            if(inter_sv_overlaps(a,b)){
                graph_put_edge(sv_graph,a,b);
                graph_put_edge(sv_graph,b,a);
            }
        }
    }
    graph_trim(sv_graph);
    vector_t *components = ic_sv_g_dfs_components(sv_graph);

    for(i=0;i<components->size;i++){
        vector_t *comp = vector_get(components,i);
        if(comp->size < MIN_INTER_CLUSTER_SIZE) { continue;}
        icclique_t *clique = icclique_find_icclique(sv_graph,comp,0,params->quasi_clique_lambda,params->quasi_clique_gamma);
        if(clique == NULL || clique->v_prime <= 0){icclique_free(clique);break;}

        vector_soft_put(calls, ic_sv_cluster_resolve(clique->items)); 
        icclique_free(clique);
    }
    vector_free(components);
    graph_free(sv_graph);
    return calls;
}

vector_t *resplit_molecules(vector_t *molecules, vector_t *discordants){
    int i;

    int base_index = -1;
    for(i=0;i<molecules->size;i++){
        interval_10X *mol = vector_get(molecules,i);
        int dis_index = interval_pair_binary_search(discordants, *mol);
        if(dis_index == -1){
            continue;
        }
        

    }
    return NULL;
}

//Direct signature -> pm, mp
//Inverted signature pp, mm
// This will be slow. 
vector_t *find_interchromosomal_events(vector_t **molecules, bam_vector_pack **reads){
	int i, j;
	sonic *snc = sonic_load(NULL);
	//parameters *params = get_params();
	vector_t *chr_to_eval = bit_set_2_index_vec( get_bam_info(NULL)->chro_bs);
	//int chr_count = params->chromosome_count;
	vector_t *all_chr_vec = vector_init(sizeof(vector_t),24*23);
    vector_t **splits = malloc(sizeof(vector_t)*get_bam_info(NULL)->chro_bs->size);
	for(j=0;j<chr_to_eval->size;j++){
		i = *(int *) vector_get(chr_to_eval,j);
		if(reads[i] == NULL){ 
			continue;
		}
		filter_dangling_reads(reads[i]->inter_pm);
		filter_dangling_reads(reads[i]->inter_mp);
		filter_dangling_reads(reads[i]->inter_pp);
		filter_dangling_reads(reads[i]->inter_mm);

        qsort(molecules[j]->items,molecules[j]->size,sizeof(void *),barcode_comp);
    }
    for(j=0;j<get_bam_info(NULL)->chro_bs->size;j++){
        if(bit_set_get_bit(get_bam_info(NULL)->chro_bs,j)){
            splits[j] = discover_split_molecules(molecules[j]);
            qsort(splits[j]->items,splits[j]->size,sizeof(void *),interval_pair_comp);
            filter_unsupported_pm_splits(splits[j],reads[j]->pm_discordants);
        }else{
            splits[j] = vector_init(sizeof(interval_pair),1);
        }
    }
	int k, t;
	for(k=0;k < chr_to_eval->size;k++){
		i = *(int *) vector_get(chr_to_eval, k);
	
/*
        for(j=0;j<reads[i]->inter_mp->size;j++){
            barcoded_read_pair *ptr = vector_get(reads[i]->inter_mp,j);
            VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\n",ptr->l_chr,ptr->left,ptr->l_or,ptr->r_chr,ptr->right,ptr->r_or,ptr->barcode);
        }
        for(j=0;j<reads[i]->inter_pm->size;j++){
            barcoded_read_pair *ptr = vector_get(reads[i]->inter_pm,j);
            VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\n",ptr->l_chr,ptr->left,ptr->l_or,ptr->r_chr,ptr->right,ptr->r_or,ptr->barcode);
        }
*/
        for(t=0;t< chr_to_eval->size;t++){
            if(k==t){continue;}
			j = *(int *) vector_get(chr_to_eval, t);
			printf("Discovering translocations between contigs %s and %s\n",snc->chromosome_names[i],snc->chromosome_names[j]);
            vector_t *pm_seps = find_inter_split_molecules(reads[i]->inter_pm,molecules[i],molecules[j]);
 //           printf("pm- %d\n",cnt);
            int kkk;

  /*          for(kkk = 0; kkk< pm_seps->size;kkk++){
                inter_split_molecule_t *ptr = vector_get(pm_seps,kkk);
                VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\t+-\n",ptr->chr1,ptr->start1,ptr->end1,ptr->chr2,ptr->start2,ptr->end2,ptr->barcode);
            }
*/

			vector_t *mp_seps = find_inter_split_molecules(reads[i]->inter_mp,molecules[i],molecules[j]);
/*			
            printf("mp- %d\n",cnt);

            for(kkk = 0; kkk< mp_seps->size;kkk++){
                inter_split_molecule_t *ptr = vector_get(mp_seps,kkk);
                VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\t-+\n",ptr->chr1,ptr->start1,ptr->end1,ptr->chr2,ptr->start2,ptr->end2,ptr->barcode);
            }
*/
						vector_t *mm_seps = find_inter_split_molecules(reads[i]->inter_mm,molecules[i],molecules[j]);
            			vector_t *pp_seps = find_inter_split_molecules(reads[i]->inter_pp,molecules[i],molecules[j]);

			vector_t *direct_t =  find_interc_translocations(pm_seps,mp_seps,splits,SV_TRANSLOCATION);
			vector_t *invert_t =  find_interc_translocations(pp_seps,mm_seps,splits,SV_INVERTED_TRANSLOCATION);

            vector_t *direct_calls = cluster_interchromosomal_events(direct_t);
            vector_t *invert_calls = cluster_interchromosomal_events(invert_t);


            vector_filter(direct_calls,ic_sv_is_proper);
            vector_filter(invert_calls,ic_sv_is_proper);
            

            vector_free(direct_t);
            vector_free(invert_t);

            vector_free(pm_seps);
            vector_free(mp_seps);
            vector_free(pp_seps);
            vector_free(mm_seps);

            vector_soft_put(all_chr_vec,direct_calls);
            vector_soft_put(all_chr_vec,invert_calls);
			//TODO
			//			vector_t *invert_t =  find_invert_translocations(ppmm_seps,mmpp_seps);
		}
	}
	return all_chr_vec;
}


