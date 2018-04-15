#include "readbam.h"
#include "common.h"
//#include "processbam.h"
#include "valorconfig.h"
#include "sonic.h"
#include "hashtable.h"

#define INITIAL_ARRAY_SIZE 500000



int altcomp(const void *v1, const void *v2){

	alt_read *r1 = *(void **)v1;
	alt_read *r2 = *(void **)v2;

	return strcmp(r1->read_name,r2->read_name);
}



void *tokenize_commas(void *val){
	return dang_string_tokenize(val,",");
}
void print_alignment_core(bam1_core_t *bam_alignment_core){
	printf( "%d\t%d\t%u\t%u\t%u\t%u\t%u\t%d\t%d\t%d\t%d\n", 
		bam_alignment_core->tid ,	
		bam_alignment_core->pos ,	
		bam_alignment_core->bin ,		
		bam_alignment_core->qual ,
		bam_alignment_core->l_qname ,
		bam_alignment_core->flag ,
		bam_alignment_core->n_cigar ,
		bam_alignment_core->l_qseq ,
		bam_alignment_core->mtid ,
		bam_alignment_core->mpos ,
		bam_alignment_core->isize);
}


bam_vector_pack *make_bam_vector_pack(){
	bam_vector_pack *new_pack = malloc(sizeof(bam_vector_pack));
	new_pack->concordants= vector_init(sizeof(interval_10X),
		INITIAL_ARRAY_SIZE);
	new_pack->mm_discordants=vector_init(sizeof(interval_discordant),
		INITIAL_ARRAY_SIZE);
	new_pack->pp_discordants=vector_init(sizeof(interval_discordant),
		INITIAL_ARRAY_SIZE);
	new_pack->pm_discordants=vector_init(sizeof(interval_discordant),
		INITIAL_ARRAY_SIZE);
;
	new_pack->mp_discordants=vector_init(sizeof(interval_discordant),
		INITIAL_ARRAY_SIZE);
;
	return new_pack;

}

bam_stats *calculate_bam_statistics( bam_info* in_bam, char* bam_path, int number_of_reads_to_stat){
	bam_stats *statistics = getMem(sizeof(bam_stats));
	
	htsFile *bam_file;

	bam_file = safe_hts_open( bam_path, "r");
	bam_hdr_t *bam_header;
	bam_header = bam_hdr_read( (bam_file->fp).bgzf);

	bam1_core_t *bam_alignment_core;
	bam1_t* bam_alignment;	
	int return_value;
	int i;
	long sum;
	vector_t *read_vec = vector_init(sizeof(int),number_of_reads_to_stat);
	
	bam_alignment = bam_init1();

	return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);

	sum = 0;
	for(i=0;i<number_of_reads_to_stat && return_value != -1 ;i++){
		bam_alignment_core = &bam_alignment->core;
	//	int flag = bam_alignment_core->flag;
		 if(    bam_alignment_core->pos != -1
			&& bam_alignment_core->mpos !=-1
			&& bam_alignment_core->tid !=-1
			&& bam_alignment_core->tid == bam_alignment_core->mtid 	
			&& bam_alignment_core->isize < 1000
			&& bam_alignment_core->isize > 0
			&& (((!bam_is_rev(bam_alignment) && bam_is_mrev(bam_alignment))||
				(bam_is_rev(bam_alignment) && !bam_is_mrev(bam_alignment)				)))){
			sum+= bam_alignment_core->isize;
			vector_put(read_vec,&(bam_alignment_core->isize));	
		}
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	}
	double mean = (double)sum/read_vec->size;
	sum = 0;
	for(i=0;i<read_vec->size;i++){
		int *val = vector_get(read_vec,i);
		sum+=(*val-mean)*(*val-mean);
	}
	double std_dev = sqrt((double)sum/read_vec->size);
	statistics->read_length_mean = mean;
	statistics->read_length_std_dev = std_dev;
	bam_hdr_destroy(bam_header);
	bam_destroy1(bam_alignment);
	if(hts_close(bam_file)){
		fprintf(stderr,"Error closing Bam file\n");
	}
	vector_free(read_vec);
	return statistics;
}

void sstrcpy ( void *dest, const void *source, size_t dummy){
	strcpy(dest,source);
}

int sstrcmp( const void *v1, const void *v2, size_t dummy){
	return strcmp(v1,v2);
}

int chr_to_tid(hashtable_t *table,char *name){
	int *tid = ht_get_value(table,name);
	return tid==NULL?-1:*tid;
}

void free_alt_read(void *vread){
	alt_read *read = vread;
	free(read->read_name);
	free(read->positions);
	free(read);
}


bam_vector_pack *read_10X_chr( bam_info* in_bam, char* bam_path, sonic *snc, int chr, bam_stats *statistics){

	int i;

//	INIT chrname Lookup table
	hashtable_t *chr_id_table = ht_init(48,sizeof(char *),sizeof(int));
	chr_id_table->key_cmp = sstrcmp;
	chr_id_table->hf = SuperFastStringHash;
	for(i=0;i<snc->number_of_chromosomes;i++){
		int *val = ht_soft_put(chr_id_table,snc->chromosome_names[i]);
		*val = i;
	}	
//


	double frag_min = MAX(0, statistics->read_length_mean - 3 * statistics->read_length_std_dev);
	double frag_max = statistics->read_length_mean + 3 * statistics->read_length_std_dev;
	static htsFile *bam_file = NULL;
	if(bam_file==NULL){
		bam_file = safe_hts_open( bam_path, "r");
	}
	static bam_hdr_t *bam_header = NULL;
	if(bam_header == NULL){
		bam_header = bam_hdr_read( (bam_file->fp).bgzf);
	}
	bam1_core_t *bam_alignment_core;

	static int return_value = 0;
	

	bam_vector_pack *pack = make_bam_vector_pack();

	static bam1_t* bam_alignment = NULL;
	

	if(bam_alignment ==NULL){
		bam_alignment = bam_init1();

	}
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);

	while(  return_value != -1 && chr == bam_alignment->core.tid){
		
		bam_alignment_core = &bam_alignment->core;
		//print_alignment_core(bam_alignment_core);
		if(bam_alignment_core->tid==-1) goto skip; //Skip invalid reads
		if(bam_alignment_core->tid >= snc->number_of_chromosomes) goto skip;
		if( bam_alignment_core->pos == -1) goto skip;
		if( bam_alignment_core->isize < 0) goto skip;
		if( bam_alignment_core->mpos == -1) goto skip;
		unsigned char * b_text =  bam_aux_get(bam_alignment,"BX");


		unsigned long barcode = encode_ten_x_barcode(b_text);

		if(barcode==-1){goto skip;}
		if(bam_alignment_core->tid!=bam_alignment_core->mtid) goto skip;

		int ccval = (is_concordant(*bam_alignment_core, frag_min, frag_max));
		int start1,end1,start2,end2;
		
		start1 = MIN(bam_alignment_core->mpos,bam_alignment_core->pos);
		start2 = MAX(bam_alignment_core->mpos,bam_alignment_core->pos);
		end1 = start1+bam_alignment_core->l_qseq;
		end2 = start2+bam_alignment_core->l_qseq;



		in_bam->read_count+=2;


		if( bam_alignment_core->qual < MIN_QUAL) goto skip;
		switch(ccval){
		case RPCONC:
			vector_put(pack->concordants,
				&(interval_10X){start1,end2,barcode}
				);
		break;
		case RPPP:
			vector_put(pack->pp_discordants,
				&(interval_discordant){start1,end1,start2,end2,barcode}
				);
		break;
		case RPMM:
			vector_put(pack->mm_discordants,
				&(interval_discordant){start1,end1,start2,end2,barcode}
				);
		break;
		case RPTDUPPM:
			vector_put(pack->pm_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
		break;
		case RPTDUPMP:
			vector_put(pack->mp_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
		break;

		default:
		break;

		}
		skip:
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	}
//	bam_destroy1(bam_alignment);
	return pack;
}



bam_vector_pack **read_10X_bam_RP( bam_info* in_bam, char* bam_path, sonic *snc){

	int i;
	bam_stats *statistics = calculate_bam_statistics(in_bam,bam_path,READ_SAMPLE_SIZE);

	vector_t *alt_reads = NULL;
	if(CHECK_ALTERNATIVE_MAPPINGS){
		vector_init(sizeof(alt_read),INITIAL_ARRAY_SIZE);
		alt_reads->rmv = free_alt_read;
	}

//	INIT chrname Lookup table
//	==================================================================
	hashtable_t *chr_id_table = ht_init(48,sizeof(char *),sizeof(int));
	chr_id_table->key_cmp = sstrcmp;
	chr_id_table->hf = SuperFastStringHash;
	for(i=0;i<snc->number_of_chromosomes;i++){
		int *val = ht_soft_put(chr_id_table,snc->chromosome_names[i]);
		*val = i;
	}	
//=========================================================================//


	double frag_min = MAX(0, statistics->read_length_mean - 3 * statistics->read_length_std_dev);
	double frag_max = statistics->read_length_mean + 3 * statistics->read_length_std_dev;
	printf("Concordant fragments are between %lf and %lf\n",frag_min, frag_max);
	htsFile *bam_file;

	bam_file = safe_hts_open( bam_path, "r");
	bam_hdr_t *bam_header;
	bam_header = bam_hdr_read( (bam_file->fp).bgzf);

	bam1_core_t *bam_alignment_core;
	bam1_t* bam_alignment;	
	int return_value;

	fprintf(stderr,"Reading BAM file %s.\n", bam_path);

	
	bam_vector_pack **packs = malloc(sizeof(bam_vector_pack)*snc->number_of_chromosomes);
	for( i = 0; i< snc->number_of_chromosomes;i++){
		packs[i] = make_bam_vector_pack();
	}
	bam_alignment = bam_init1();
	return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);

	while(  return_value != -1){
		
		bam_alignment_core = &bam_alignment->core;
		//print_alignment_core(bam_alignment_core);
		if(bam_alignment_core->tid==-1) goto skip; //Skip invalid reads
		if(bam_alignment_core->tid >= snc->number_of_chromosomes) goto skip;
		if( bam_alignment_core->pos == -1) goto skip;

		if( bam_alignment_core->mpos == -1) goto skip;
		unsigned long barcode = encode_ten_x_barcode( bam_aux_get(bam_alignment,"BX"));
		if( CHECK_ALTERNATIVE_MAPPINGS){

				alt_read *new_alt = malloc(sizeof(alt_read));
				char * read_name = bam_get_qname(bam_alignment);
				new_alt->read_name = malloc(sizeof(char) * (strlen(read_name) + 1));
				strcpy(new_alt->read_name,read_name);
				new_alt->barcode = barcode;

				if(barcode==-1){goto skip;}
				unsigned char *flag = bam_aux_get(bam_alignment,ALTERNATIVE_MAPPING_FLAG);
				if(flag!=NULL){
					flag+=1;//Skip the Z
				}
		// ALTERNATIVE MAPPING string should tokenize to 6 pieces, otherwise
		// crash
				vector_t *alts = dang_string_tokenize((char *)flag,";,");

				new_alt->count = 1+alts->size/6;

				new_alt->positions = malloc(sizeof(simple_interval) * new_alt->count);
				new_alt->flag =  bam_alignment_core->flag;
				new_alt->positions[0].tid = bam_alignment_core->tid;
				new_alt->positions[0].start = bam_alignment_core->pos;
				new_alt->positions[0].end = bam_alignment_core->pos + bam_alignment_core->l_qseq;
				new_alt->positions[0].strand = bam_is_rev(bam_alignment); 

				int new_count = new_alt->count;
				alts->REMOVE_POLICY= REMP_LAZY;
				for(i=0;i<new_alt->count-1;i++){
					int qual = atoi(vector_get(alts,6*i+4));
					int tid = chr_to_tid(chr_id_table,vector_get(alts,6*i));
					if( tid == -1 || 0.75 * bam_alignment_core->qual  > qual+1){
						new_count--;
						int k;
						for(k=0;k<6;k++){
							vector_remove(alts,6*i+k);
						}
					}
				}
				vector_defragment(alts);
				alts->REMOVE_POLICY = REMP_SORTED;
				new_alt->count = new_count;
				for(i=0;i<new_alt->count-1;i++){
					new_alt->positions[1+i].tid = chr_to_tid(chr_id_table,vector_get(alts,6*i));
					new_alt->positions[1+i].start = atoi(vector_get(alts,6*i+1));
					new_alt->positions[1+i].end = bam_alignment_core->l_qseq + new_alt->positions[1+i].start;
					new_alt->positions[1+i].strand = strcmp(vector_get(alts,6*i+2),"+")==0?READ_STRAND_POS:READ_STRAND_NEG;
				}
				
				vector_free(alts);
				vector_soft_put(alt_reads,new_alt);
			
		}else{

			int ccval = (is_concordant(*bam_alignment_core, frag_min, frag_max));
			int start1,end1,start2,end2;
			
			start1 = MIN(bam_alignment_core->mpos,bam_alignment_core->pos);
			start2 = MAX(bam_alignment_core->mpos,bam_alignment_core->pos);
			end1 = start1+bam_alignment_core->l_qseq;
			end2 = start2+bam_alignment_core->l_qseq;



			in_bam->read_count+=2;

			if(barcode==-1){goto skip;}
			if( bam_alignment_core->qual < MIN_QUAL) goto skip;
			switch(ccval){
			case RPCONC:
				vector_put(packs[bam_alignment_core->tid]->concordants,
					&(interval_10X){start1,end2,barcode}
					);
			break;
			case RPPP:
				vector_put(packs[bam_alignment_core->tid]->pp_discordants,
					&(interval_discordant){start1,end1,start2,end2,barcode}
					);
			break;
			case RPMM:
				vector_put(packs[bam_alignment_core->tid]->mm_discordants,
					&(interval_discordant){start1,end1,start2,end2,barcode}
					);
			break;
			case RPTDUPPM:
				vector_put(packs[bam_alignment_core->tid]->pm_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
			break;
			case RPTDUPMP:
				vector_put(packs[bam_alignment_core->tid]->mp_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
			break;

			default:
			break;
			}
		}
		skip:
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	}
	if(CHECK_ALTERNATIVE_MAPPINGS){
//This might take some time...
		qsort(alt_reads->items,alt_reads->size,sizeof( void *),altcomp);
		for(i=1;i<alt_reads->size;i++){
			alt_read *a1 = vector_get(alt_reads,i-1);
			alt_read *a2 = vector_get(alt_reads,i);
			if(strcmp(a1->read_name,a2->read_name)==0){
				int j,k;
				for(j=0;j<a1->count;j++){
					for(k=0;k<a2->count;k++){
						if(a1->positions[j].tid != a2->positions[k].tid) continue;

						int ccval = (is_alt_concordant(a1->positions[j].start,a2->positions[k].end,a1->flag,
								a2->positions[k].strand, a1->positions[j].strand, frag_min, frag_max));
						int start1,end1,start2,end2;
						

						if( a1->positions[j].start > a2->positions[k].start){
							start1 = a2->positions[k].start;
							start2 = a1->positions[j].start;;
							end1 = a2->positions[k].end;
							end2 = a1->positions[j].end;
						}
						else{
							start2 = a2->positions[k].start;
							start1 = a1->positions[j].start;;
							end2 = a2->positions[k].end;
							end1 = a1->positions[j].end;
						}
						switch(ccval){
						case RPCONC:
							vector_put(packs[a1->positions[j].tid]->concordants,
								&(interval_10X){start1,end2,a1->barcode}
								);
						break;
						case RPPP:
							vector_put(packs[a1->positions[j].tid]->pp_discordants,
								&(interval_discordant){start1,end1,start2,end2,a1->barcode}
								);
						break;
						case RPMM:
							vector_put(packs[a1->positions[j].tid]->mm_discordants,
								&(interval_discordant){start1,end1,start2,end2,a1->barcode}
								);
						break;
						case RPTDUPPM:
							vector_put(packs[a1->positions[j].tid]->pm_discordants,&(interval_discordant){start1,end1,start2,end2,a1->barcode});
						break;
						case RPTDUPMP:
							vector_put(packs[a1->positions[j].tid]->mp_discordants,&(interval_discordant){start1,end1,start2,end2,a1->barcode});
						break;

						default:
						break;
						}
				
					}
				}
				i++;
			}
		}
	}
	free(statistics);
	vector_free(alt_reads);
	if(hts_close(bam_file)){
		fprintf(stderr,"Error closing Bam file\n");
	}
	
	bam_hdr_destroy(bam_header);
	bam_destroy1(bam_alignment);
	return packs;
}



bam_vector_pack **read_10X_bam( bam_info* in_bam, char* bam_path, sonic *snc){

	int i;
	bam_stats *statistics = calculate_bam_statistics(in_bam,bam_path,READ_SAMPLE_SIZE);

	vector_t *alt_reads = vector_init(sizeof(alt_read),INITIAL_ARRAY_SIZE);
	alt_reads->rmv = free_alt_read;
//	INIT chrname Lookup table
	hashtable_t *chr_id_table = ht_init(48,sizeof(char *),sizeof(int));
	chr_id_table->key_cmp = sstrcmp;
	chr_id_table->hf = SuperFastStringHash;
	for(i=0;i<snc->number_of_chromosomes;i++){
		int *val = ht_soft_put(chr_id_table,snc->chromosome_names[i]);
		*val = i;
	}	
//


	double frag_min = MAX(0, statistics->read_length_mean - 3 * statistics->read_length_std_dev);
	double frag_max = statistics->read_length_mean + 3 * statistics->read_length_std_dev;
	printf("Concordant fragments are between %lf and %lf\n",frag_min, frag_max);
	htsFile *bam_file;

	bam_file = safe_hts_open( bam_path, "r");
	bam_hdr_t *bam_header;
	bam_header = bam_hdr_read( (bam_file->fp).bgzf);

	bam1_core_t *bam_alignment_core;
	bam1_t* bam_alignment;	
	int return_value;

	fprintf(stderr,"Reading BAM file %s.\n", bam_path);

	
	bam_vector_pack **packs = malloc(sizeof(bam_vector_pack)*snc->number_of_chromosomes);
	for( i = 0; i< snc->number_of_chromosomes;i++){
		packs[i] = make_bam_vector_pack();
	}
	bam_alignment = bam_init1();
	return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);

	while(  return_value != -1){
		
		bam_alignment_core = &bam_alignment->core;
		//print_alignment_core(bam_alignment_core);
		if(bam_alignment_core->tid==-1) goto skip; //Skip invalid reads
		if(bam_alignment_core->tid >= snc->number_of_chromosomes) goto skip;
		if( bam_alignment_core->pos == -1) goto skip;

		if( bam_alignment_core->mpos == -1) goto skip;
		unsigned long barcode = encode_ten_x_barcode( bam_aux_get(bam_alignment,"BX"));
		if( CHECK_ALTERNATIVE_MAPPINGS){

				alt_read *new_alt = malloc(sizeof(alt_read));
				char * read_name = bam_get_qname(bam_alignment);
				new_alt->read_name = malloc(sizeof(char) * (strlen(read_name) + 1));
				strcpy(new_alt->read_name,read_name);
				new_alt->barcode = barcode;

				if(barcode==-1){goto skip;}
				unsigned char *flag = bam_aux_get(bam_alignment,ALTERNATIVE_MAPPING_FLAG);
				if(flag!=NULL){
					flag+=1;//Skip the Z
				}
				vector_t *alts = dang_string_tokenize((char *)flag,";,");

				new_alt->count = 1+alts->size/6;

				new_alt->positions = malloc(sizeof(simple_interval) * new_alt->count);
				new_alt->flag =  bam_alignment_core->flag;
				new_alt->positions[0].tid = bam_alignment_core->tid;
				new_alt->positions[0].start = bam_alignment_core->pos;
				new_alt->positions[0].end = bam_alignment_core->pos + bam_alignment_core->l_qseq;
				new_alt->positions[0].strand = bam_is_rev(bam_alignment); 

				int new_count = new_alt->count;
				alts->REMOVE_POLICY= REMP_LAZY;
				for(i=0;i<new_alt->count-1;i++){
					int qual = atoi(vector_get(alts,6*i+4));
					int tid = chr_to_tid(chr_id_table,vector_get(alts,6*i));
					if( tid == -1 || 0.75 * bam_alignment_core->qual  > qual+1){
						new_count--;
						int k;
						for(k=0;k<6;k++){
							vector_remove(alts,6*i+k);
						}
					}
				}
				vector_defragment(alts);
				alts->REMOVE_POLICY = REMP_SORTED;
				new_alt->count = new_count;
				for(i=0;i<new_alt->count-1;i++){
					new_alt->positions[1+i].tid = chr_to_tid(chr_id_table,vector_get(alts,6*i));
					new_alt->positions[1+i].start = atoi(vector_get(alts,6*i+1));
					new_alt->positions[1+i].end = bam_alignment_core->l_qseq + new_alt->positions[1+i].start;
					new_alt->positions[1+i].strand = strcmp(vector_get(alts,6*i+2),"+")==0?READ_STRAND_POS:READ_STRAND_NEG;
				}
				
				vector_free(alts);
				vector_soft_put(alt_reads,new_alt);
			
		}else{

			if(bam_alignment_core->tid!=bam_alignment_core->mtid) goto skip;

			int ccval = (is_concordant(*bam_alignment_core, frag_min, frag_max));
			int start1,end1,start2,end2;
			
			start1 = MIN(bam_alignment_core->mpos,bam_alignment_core->pos);
			start2 = MAX(bam_alignment_core->mpos,bam_alignment_core->pos);
			end1 = start1+bam_alignment_core->l_qseq;
			end2 = start2+bam_alignment_core->l_qseq;



			in_bam->read_count+=2;

			if(barcode==-1){goto skip;}
			if( bam_alignment_core->qual < MIN_QUAL) goto skip;
			switch(ccval){
			case RPCONC:
				vector_put(packs[bam_alignment_core->tid]->concordants,
					&(interval_10X){start1,end2,barcode}
					);
			break;
			case RPPP:
				vector_put(packs[bam_alignment_core->tid]->pp_discordants,
					&(interval_discordant){start1,end1,start2,end2,barcode}
					);
			break;
			case RPMM:
				vector_put(packs[bam_alignment_core->tid]->mm_discordants,
					&(interval_discordant){start1,end1,start2,end2,barcode}
					);
			break;
			case RPTDUPPM:
				vector_put(packs[bam_alignment_core->tid]->pm_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
			break;
			case RPTDUPMP:
				vector_put(packs[bam_alignment_core->tid]->mp_discordants,&(interval_discordant){start1,end1,start2,end2,barcode});
			break;

			default:
			break;
			}
		}
		skip:
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	}
	if(CHECK_ALTERNATIVE_MAPPINGS){
//This might take some time...
		qsort(alt_reads->items,alt_reads->size,sizeof( void *),altcomp);
		for(i=1;i<alt_reads->size;i++){
			alt_read *a1 = vector_get(alt_reads,i-1);
			alt_read *a2 = vector_get(alt_reads,i);
			if(strcmp(a1->read_name,a2->read_name)==0){
				int j,k;
				for(j=0;j<a1->count;j++){
					for(k=0;k<a2->count;k++){
						if(a1->positions[j].tid != a2->positions[k].tid) continue;

						int ccval = (is_alt_concordant(a1->positions[j].start,a2->positions[k].end,a1->flag,
								a2->positions[k].strand, a1->positions[j].strand, frag_min, frag_max));
						int start1,end1,start2,end2;
						

						if( a1->positions[j].start > a2->positions[k].start){
							start1 = a2->positions[k].start;
							start2 = a1->positions[j].start;;
							end1 = a2->positions[k].end;
							end2 = a1->positions[j].end;
						}
						else{
							start2 = a2->positions[k].start;
							start1 = a1->positions[j].start;;
							end2 = a2->positions[k].end;
							end1 = a1->positions[j].end;
						}
						switch(ccval){
						case RPCONC:
							vector_put(packs[a1->positions[j].tid]->concordants,
								&(interval_10X){start1,end2,a1->barcode}
								);
						break;
						case RPPP:
							vector_put(packs[a1->positions[j].tid]->pp_discordants,
								&(interval_discordant){start1,end1,start2,end2,a1->barcode}
								);
						break;
						case RPMM:
							vector_put(packs[a1->positions[j].tid]->mm_discordants,
								&(interval_discordant){start1,end1,start2,end2,a1->barcode}
								);
						break;
						case RPTDUPPM:
							vector_put(packs[a1->positions[j].tid]->pm_discordants,&(interval_discordant){start1,end1,start2,end2,a1->barcode});
						break;
						case RPTDUPMP:
							vector_put(packs[a1->positions[j].tid]->mp_discordants,&(interval_discordant){start1,end1,start2,end2,a1->barcode});
						break;

						default:
						break;
						}
				
					}
				}
				i++;
			}
		}
	}
	free(statistics);
	vector_free(alt_reads);
	if(hts_close(bam_file)){
		fprintf(stderr,"Error closing Bam file\n");
	}
	
	bam_hdr_destroy(bam_header);
	bam_destroy1(bam_alignment);
	return packs;
}


void destroy_bams(bam_vector_pack *reads){
	vector_free(reads->concordants);
	vector_free(reads->pp_discordants);
	vector_free(reads->mm_discordants);
	vector_free(reads->pm_discordants);
	vector_free(reads->mp_discordants);
	free(reads);
}
