#include "readbam.h"
#include "common.h"
//#include "processbam.h"
#include "valorconfig.h"
#include "sonic.h"
#define INITIAL_ARRAY_SIZE 25000
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


bam_vector_pack *make_bam_vector_pack(int count){
	int i;
	vector_t **vectors = (vector_t **) getMem(sizeof(vector_t *) * count); //A vector of intervals for each alignment	
	vector_t **ppvectors = (vector_t **) getMem(sizeof(vector_t *) * count); //A vector of intervals for each alignment
	vector_t **mmvectors = (vector_t **) getMem(sizeof(vector_t *) * count); //A vector of intervals for each alignment
	vector_t **pmdup = (vector_t **) getMem(sizeof(vector_t *) * count); //A vector of intervals for each alignment
	vector_t **mpdup = (vector_t **) getMem(sizeof(vector_t *) * count); //A vector of intervals for each alignment
	
	for( i=0; i<count;i++){
		vectors[i] = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
		mmvectors[i] = vector_init(sizeof(interval_discordant),INITIAL_ARRAY_SIZE);
		ppvectors[i] = vector_init(sizeof(interval_discordant),INITIAL_ARRAY_SIZE);
		pmdup[i] = vector_init(sizeof(interval_discordant),INITIAL_ARRAY_SIZE);
		mpdup[i] = vector_init(sizeof(interval_discordant),INITIAL_ARRAY_SIZE);
		
	}
	bam_vector_pack *new_pack = malloc(sizeof(bam_vector_pack));
	new_pack->concordants=vectors;
	new_pack->mm_discordants=mmvectors;
	new_pack->pp_discordants=ppvectors;
	new_pack->pm_dup=pmdup;
	new_pack->mp_dup=mpdup;
	new_pack->count = count;
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
//	vector_free(read_vec);
	return statistics;
}


bam_vector_pack *read_10X_bam( bam_info* in_bam, char* bam_path, sonic *snc){

	int i;
	bam_stats *statistics = calculate_bam_statistics(in_bam,bam_path,READ_SAMPLE_SIZE);


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
	in_bam->read_depth = getMem(sizeof(short*)*bam_header->n_targets);
	in_bam->chromosome_mean_rd = getMem(sizeof(double)*bam_header->n_targets);
	for(i=0;i<snc->number_of_chromosomes;i++){
		in_bam->read_depth[i]=getMem(sizeof(short)*snc->chromosome_lengths[i]);
		int j;
		for( j=0;j<snc->chromosome_lengths[i];j++){
			in_bam->read_depth[i][j]=0;
		}
	}	
	
	fprintf(stderr,"Reading BAM file %s.\n", bam_path);

	
	bam_vector_pack *pack = make_bam_vector_pack(bam_header->n_targets);

	bam_alignment = bam_init1();
	return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	int valid = 0;
	int counter = 0;
	int ppcnt = 0;
	int mmcnt = 0;
	while(  return_value != -1){
		
		bam_alignment_core = &bam_alignment->core;
		//print_alignment_core(bam_alignment_core);
		if(bam_alignment_core->tid==-1) goto skip; //Skip invalid reads
		if(bam_alignment_core->tid!=bam_alignment_core->mtid) goto skip;
//		if( bam_alignment_core->isize < 0) goto skip;

		if( bam_alignment_core->pos == -1) goto skip;


		int ccval = (is_concordant(*bam_alignment_core, frag_min, frag_max));
		int start1,end1,start2,end2;

		unsigned long barcode = encode_ten_x_barcode( bam_aux_get(bam_alignment,"BX"));

//		printf("Barcode: %lu\n",barcode);
		start1 =MIN(bam_alignment_core->mpos,bam_alignment_core->pos);
		start2 =MAX(bam_alignment_core->mpos,bam_alignment_core->pos);
		end1=start1+bam_alignment_core->l_qseq;
		end2=start2+bam_alignment_core->l_qseq;


		in_bam->read_depth[bam_alignment_core->tid][
			bam_alignment_core->pos]++;

		if( bam_alignment_core->mpos == -1) goto skip;
		in_bam->read_depth[bam_alignment_core->tid][
			bam_alignment_core->mpos]++;
		in_bam->read_count+=2;


		if(barcode==-1){goto skip;}
		if( bam_alignment_core->qual < MIN_QUAL) goto skip;
		switch(ccval){
		case RPCONC:
		//			in_bam->read_depth[bam_alignment_core->tid][end2]++;
			valid ++;
			vector_put(pack->concordants[bam_alignment_core->tid],
				&(interval_10X){start1,end2,barcode}
				);
		break;
		case RPPP:
			vector_put(pack->pp_discordants[bam_alignment_core->tid],
				&(interval_discordant){start1,end1,start2,end2,barcode}
				);
			ppcnt ++;		
		break;
		case RPMM:
			vector_put(pack->mm_discordants[bam_alignment_core->tid],
				&(interval_discordant){start1,end1,start2,end2,barcode}
				);
			mmcnt ++;
		break;
		case RPTDUPPM:

			VALOR_LOG("%d %d %d %d +-\n",start1,end1,start2,end2);
			vector_put(pack->pm_dup[bam_alignment_core->tid],&(interval_discordant){start1,end1,start2,end2,barcode});
		break;
		case RPTDUPMP:

			VALOR_LOG("%d %d %d %d -+\n",start1,end1,start2,end2);
			vector_put(pack->mp_dup[bam_alignment_core->tid],&(interval_discordant){start1,end1,start2,end2,barcode});
		break;

		default:
		break;
		}
		skip:
		counter ++;
		return_value = bam_read1( (bam_file->fp).bgzf, bam_alignment);
	}
	free(statistics);
	bam_hdr_destroy(bam_header);
	bam_destroy1(bam_alignment);
	if(hts_close(bam_file)){
		fprintf(stderr,"Error closing Bam file\n");
	}
	printf("Read %d valid lines out of %d reads\n++discordants:%d\t--discordants:%d\n",valid,counter,ppcnt,mmcnt);
	
	return pack;
}


void destroy_bams(bam_vector_pack *reads){
	int i;
	for( i = 0; i < reads->count ; i++){
		vector_free(reads->concordants[i]);
		vector_free(reads->pp_discordants[i]);
		vector_free(reads->mm_discordants[i]);
		vector_free(reads->pm_dup[i]);
		vector_free(reads->mp_dup[i]);
	}

	freeMem(reads->concordants,sizeof(vector_t **));
	freeMem(reads->mm_discordants,sizeof(vector_t **));
	freeMem(reads->pp_discordants,sizeof(vector_t **));
	freeMem(reads->pm_dup,sizeof(vector_t **));
	freeMem(reads->mp_dup,sizeof(vector_t **));
	freeMem(reads,sizeof(bam_vector_pack));

}
