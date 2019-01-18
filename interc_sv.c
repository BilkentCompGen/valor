#include "interc_sv.h"


#define SCL_INIT_LIMIT 400




size_t split_molecule_binary_search(vector_t *splits, interval_10X key){
	if(splits->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = splits->size - 1;

	while( first < last){
		mid = (first + last)/2;
		if(((splitmolecule_t*)vector_get(splits,mid))->start1< key.start){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
	}

	while( mid >= 0 && key.start -30000< ((splitmolecule_t*) vector_get(splits,mid))->start1){
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
size_t barcode_binary_search(vector_t *mols, unsigned long key){
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
	return mid +1;
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
int inter_split_overlaps(inter_split_molecule_t s1, inter_split_molecule_t s2){
	if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}
	return 0;
	//TODO
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


int inter_split_indicates_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2){
	if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}
    if( s1.barcode == s2.barcode) {return 0;}

	interval_10X A = (interval_10X){s1.start1,s1.end1,s1.barcode};  
	interval_10X B = (interval_10X){s1.start2,s1.end2,s1.barcode};  
	interval_10X C = (interval_10X){s2.start1,s2.end1,s2.barcode};  
	interval_10X D = (interval_10X){s2.start2,s2.end2,s2.barcode};  
	if(interval_inner_distance(A,C) < DUP_GAP &&
			interval_inner_distance(A,C) > DUP_OVERLAP){
		return INTERC_BACK_COPY;
	}
	else if(interval_inner_distance(B,D) < DUP_GAP &&
			interval_inner_distance(B,D) > DUP_OVERLAP){
		return INTERC_FORW_COPY;
	}
	return 0; 
}

ic_sv_t *inter_sv_init(inter_split_molecule_t *a, inter_split_molecule_t *b, splitmolecule_t *tra_del, sv_type type, int orient){
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

splitmolecule_t *find_matching_split_molecule(vector_t **splits,inter_split_molecule_t *a, inter_split_molecule_t *b, int orient){    
    
	if(orient == 0){
       
        VALOR_LOG("BAD Orient: %d\n",orient);
        return NULL;}
    int src_chr; 
    int tar_chr;

    interval_10X deletion_interval = {0,0,0};
    if( orient == INTERC_FORW_COPY){
        src_chr = MIN(a->chr1,a->chr2);
        tar_chr = MAX(a->chr1,a->chr2);
        deletion_interval = (interval_10X){.start = a->end1,.end = b->start1,.barcode=0};
    }else if(orient == INTERC_BACK_COPY){
        tar_chr = MIN(a->chr1,a->chr2);
        src_chr = MAX(a->chr1,a->chr2);
        deletion_interval = (interval_10X){.start = a->end2,.end = b->start2,.barcode=0};

    }else{
        
        fprintf(stderr, "Translocation orientation cannot be %d, terminating!!\n",orient);
        exit(-1);
        return NULL;
    }


    size_t pos = split_molecule_binary_search(splits[src_chr],deletion_interval);
    if( pos == -1){
        VALOR_LOG("NO %d %d\n",deletion_interval.start,deletion_interval.end);
        return NULL;
    }

    interval_pair *cand = vector_get(splits[src_chr],pos);
    while( pos < splits[src_chr]->size && cand->start1 < deletion_interval.end + 50000){

        if( cand->start1 + cand->end1 > 2 * deletion_interval.start - 50000 && 
                cand->start2 + cand->end2 < 2 * deletion_interval.end + 50000){
    
            return vector_get(splits[src_chr], pos);
        }
        pos++;
        cand = vector_get(splits[src_chr],pos);
    }

    VALOR_LOG("CM\n%d\t%d\t%d\n",src_chr,deletion_interval.start,deletion_interval.end);
    return NULL;
}

vector_t *find_direct_translocations(vector_t *sp1, vector_t *sp2, vector_t **molecules){
	vector_t *tlocs = vector_init(sizeof(ic_sv_t),512);
	int i,j;

	for(i=0;i<sp1->size;i++){
		inter_split_molecule_t *a = vector_get(sp1,i);
		for(j=0;j<sp2->size;j++){   
			inter_split_molecule_t *b = vector_get(sp2,j);
			int orient = inter_split_indicates_translocation(*a,*b);

            splitmolecule_t *del_tra = find_matching_split_molecule(molecules,a,b,orient);
			
            VALOR_LOG("%d\t%d\n",orient,del_tra!=NULL);
            if(orient && del_tra != NULL){
				vector_soft_put(tlocs,inter_sv_init(a,b,del_tra,SV_TRANSLOCATION,orient));
			}
		}
	}

	return tlocs;
}

int is_inter_chr_split(barcoded_read_pair *pair, interval_10X *a, interval_10X *b){
	if( (pair->barcode !=a->barcode) || (pair->barcode!=b->barcode)){ return 0;}

	return (in_range(pair->left,a->start,2*MOLECULE_EXT) && in_range(pair->right,b->start,2*MOLECULE_EXT));        
}

int cnt = 0;
// returns a vector of inter_split_molecules
vector_t *find_separated_molecules(vector_t *reads, vector_t *mol_a, vector_t *mol_b){
	vector_t *isms = vector_init(sizeof(inter_split_molecule_t),512);
	int i;

	for(i=0;i<reads->size;i++){
		barcoded_read_pair *pair =  vector_get(reads,i);
		if( pair->barcode == -1){ continue;}
        size_t k = 0;
		//size_t k = barcode_binary_search(mol_a,pair->barcode);
         while (k < mol_a->size && ((interval_10X *)vector_get(mol_a,k))->barcode != pair->barcode){
            k++;
        }
        if ( k >= mol_a->size){
            continue;
        }
		//size_t t = barcode_binary_search(mol_b,pair->barcode);
		size_t t = 0;
        while (t < mol_b->size && ((interval_10X *)vector_get(mol_b,t))->barcode != pair->barcode){
            t++;
        }
        if ( t >= mol_b->size){
            continue;
        }
        size_t kk = k;
		size_t tt = t;

		while(I10X_VECTOR_GET(mol_a,kk)->barcode == pair->barcode){	
			while(I10X_VECTOR_GET(mol_b,tt)->barcode == pair->barcode){

                        
                if(is_inter_chr_split(pair,vector_get(mol_a,kk),vector_get(mol_b,tt))){
					cnt++;
                    vector_soft_put(isms,inter_split_init(pair,vector_get(mol_a,kk),vector_get(mol_b,tt)));
			    }else  if(is_inter_chr_split(pair,vector_get(mol_b,tt),vector_get(mol_a,kk))){
					cnt++;
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

//Direct signature -> pm, mp
//Inverted signature pp, mm
// This will be slow. 
vector_t *find_interchromosomal_events(vector_t **molecules, bam_vector_pack **reads){
	int i, j;
	//sonic *snc = sonic_load(NULL);
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
        for(t=k+1;t< chr_to_eval->size;t++){
			j = *(int *) vector_get(chr_to_eval, t);
			vector_t *pm_seps = find_separated_molecules(reads[i]->inter_pm,molecules[i],molecules[j]);
 //           printf("pm- %d\n",cnt);
            int kkk;

  /*          for(kkk = 0; kkk< pm_seps->size;kkk++){
                inter_split_molecule_t *ptr = vector_get(pm_seps,kkk);
                VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\t+-\n",ptr->chr1,ptr->start1,ptr->end1,ptr->chr2,ptr->start2,ptr->end2,ptr->barcode);
            }
*/

			vector_t *mp_seps = find_separated_molecules(reads[i]->inter_mp,molecules[i],molecules[j]);
/*			
            printf("mp- %d\n",cnt);
            //			vector_t *pp_seps = find_separated_molecules(reads[i]->inter_pp,molecules[i],molecules[j]);
            for(kkk = 0; kkk< mp_seps->size;kkk++){
                inter_split_molecule_t *ptr = vector_get(mp_seps,kkk);
                VALOR_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%lu\t-+\n",ptr->chr1,ptr->start1,ptr->end1,ptr->chr2,ptr->start2,ptr->end2,ptr->barcode);
            }
*/
			//			vector_t *mm_seps = find_separated_molecules(reads[i]->inter_mm,molecules[i],molecules[j]);

			vector_t *direct_t =  find_direct_translocations(pm_seps,mp_seps,splits);
			vector_soft_put(all_chr_vec,direct_t);
			//TODO
			//			vector_t *invert_t =  find_invert_translocations(ppmm_seps,mmpp_seps);
		}
	}
	return all_chr_vec;
}


