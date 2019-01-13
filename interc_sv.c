#include "interc_sv.h"


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
		mid = (first + last)/2;
		if(((interval_10X *)vector_get(mols,mid))->barcode < key){
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
	fprintf(stream,"%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
			snc->chromosome_names[sv->chr_source],
			sv->AB.start1,
			sv->AB.end1,
			snc->chromosome_names[sv->chr_target],
			sv->CD.start1,
			sv->CD.start1,
			sv_type_name(sv->type)
	       );
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

	return isplit;
}


int inter_split_indicates_translocation(inter_split_molecule_t s1, inter_split_molecule_t s2){
	if( !(s1.chr1 == s2.chr1 && s1.chr2 == s2.chr2)){ return 0;}

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

ic_sv_t *inter_sv_init(inter_split_molecule_t *s1, inter_split_molecule_t *s2, sv_type type, int orient){
	ic_sv_t *new_i = malloc(sizeof(ic_sv_t));
	memset(new_i,0,sizeof(ic_sv_t));
	new_i->supports[0] = 1;//TODO fix this
	new_i->supports[1] = 1;//TODO fix this
	new_i->supports[2] = 1;//TODO fix this

	new_i->covered = 0;
	new_i->tabu = 0;
	new_i->dv = 0;
	new_i->inactive = 0;
	new_i->AB=*s1;
	new_i->CD=*s2;
	new_i->EF=(splitmolecule_t){0,0,0,0,0L};
	new_i->type=type;
	//TODO assign chr's accordign to orientation.
	new_i->chr_source = s1->chr1;
	new_i->chr_target = s1->chr2;

	return new_i;
}
vector_t *find_direct_translocations(vector_t *sp1, vector_t *sp2){
	vector_t *tlocs = vector_init(sizeof(ic_sv_t),512);
	int i,j;

	for(i=0;i<sp1->size;i++){
		inter_split_molecule_t *a = vector_get(sp1,i);
		for(j=0;j<sp2->size;j++){   
			inter_split_molecule_t *b = vector_get(sp2,j);
			int orient = inter_split_indicates_translocation(*a,*b);
			if(orient){
				vector_soft_put(tlocs,inter_sv_init(a,b,SV_TRANSLOCATION,orient));
			}
		}
	}

	return tlocs;
}

int is_inter_chr_split(barcoded_read_pair *pair, interval_10X *a, interval_10X *b){
	if( (pair->barcode !=a->barcode) || (pair->barcode!=b->barcode)){ return 0;}
	return in_range(pair->left,a->start,2*MOLECULE_EXT) && in_range(pair->right,b->start,2*MOLECULE_EXT);
}

// returns a vector of inter_split_molecules
vector_t *find_separated_molecules(vector_t *reads, vector_t *mol_a, vector_t *mol_b){
	vector_t *isms = vector_init(sizeof(inter_split_molecule_t),512);
	int i;

	for(i=0;i<reads->size;i++){
		barcoded_read_pair *pair =  vector_get(reads,i);
		if( pair->barcode == -1){ continue;}
		size_t k = barcode_binary_search(mol_a,pair->barcode);
		size_t t = barcode_binary_search(mol_b,pair->barcode);
		size_t kk = k;
		size_t tt = t;
		while(I10X_VECTOR_GET(mol_a,kk)->barcode == pair->barcode){	
			while(I10X_VECTOR_GET(mol_b,tt)->barcode == pair->barcode){
				if(is_inter_chr_split(pair,vector_get(mol_a,kk),vector_get(mol_b,tt))){
					vector_soft_put(isms,inter_split_init(pair,vector_get(mol_a,kk),vector_get(mol_b,tt)));
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

	for(j=0;j<chr_to_eval->size;j++){
		i = *(int *) vector_get(chr_to_eval,j);
		if(reads[i] == NULL){ 
			continue;
		}
		filter_dangling_reads(reads[i]->inter_pm);
		filter_dangling_reads(reads[i]->inter_mp);
		filter_dangling_reads(reads[i]->inter_pp);
		filter_dangling_reads(reads[i]->inter_mm);

	}
	int k, t;
	for(k=0;k < chr_to_eval->size;k++){
		i = *(int *) vector_get(chr_to_eval, k);
		for(t=k+1;t< chr_to_eval->size;t++){
			j = *(int *) vector_get(chr_to_eval, t);
			vector_t *pm_seps = find_separated_molecules(reads[i]->inter_pm,molecules[i],molecules[j]);
			vector_t *mp_seps = find_separated_molecules(reads[i]->inter_mp,molecules[i],molecules[j]);
			//			vector_t *pp_seps = find_separated_molecules(reads[i]->inter_pp,molecules[i],molecules[j]);
			//			vector_t *mm_seps = find_separated_molecules(reads[i]->inter_mm,molecules[i],molecules[j]);

			vector_t *direct_t =  find_direct_translocations(pm_seps,mp_seps);
			vector_soft_put(all_chr_vec,direct_t);
			//TODO
			//			vector_t *invert_t =  find_invert_translocations(ppmm_seps,mmpp_seps);
		}
	}
	return all_chr_vec;
}


