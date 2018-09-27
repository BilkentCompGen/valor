#include "cnv.h"



double get_depth_region(short *depths, int start, int end){
	if( end <start){ return get_depth_region(depths,end,start);}
	double sum = 0;
	int i;
	
	for( i = floor(start/MOLECULE_BIN_SIZE);i<ceil(end/MOLECULE_BIN_SIZE);i++){
		sum+=depths[i];
	}
	return sum / ceil((1.0*end-start)/MOLECULE_BIN_SIZE);
}

double get_depth_deviation(short *depths, int start, int end){
	if( end <start){ return get_depth_deviation(depths,end,start);}
    double mean = get_depth_region(depths,start,end);

    int i;
    double sum = 0;
    
    for(i= floor(start/MOLECULE_BIN_SIZE);i<ceil(end/MOLECULE_BIN_SIZE);i++){
        sum+=((depths[i]-mean) * (depths[i]-mean));
    }
    return sqrt(sum/ceil((1.0*end-start)/MOLECULE_BIN_SIZE));
}

int cmp_short(const void *a, const void *b){
	short aa = *(short *) a;
	short bb = *(short *) b;
	return aa - bb;
}

/*
double calculate_chromosome_RD_mean(bam_info *in_bam, sonic *snc, int chr_no){

	double mean;
	int i;
	long sum;
	
	int chr_len = 	snc->chromosome_lengths[chr_no];
	
	short *read_depths = getMem(sizeof(short)*chr_len);
 	memcpy(read_depths, in_bam->read_depth[chr_no],sizeof(short)*chr_len);
	

	
	sum = 0;
	for( i=0;i<chr_len;i++){
		sum+=read_depths[i];	
	}
	mean = (double)sum/(chr_len);
	printf("Initial Chromosome %s mean is %lf\n",snc->chromosome_names[chr_no],mean);

 // Followinw Code should remove the outlier read depths
 // not worth for the time it uses
//
	qsort(read_depths,chr_len,sizeof(short),cmp_short);
	sum = 0;
	for( i=0;i<chr_len;i++){
		sum+=(read_depths[i]-mean)*(read_depths[i]-mean);
	}
	std_dev = sqrt((double)sum/chr_len);
	while( std_dev > ( mean/4)){
	
//		if( abs(read_depths[start]-mean) > abs(read_depths[chr_len-1]-mean)){
//			mean= (mean*(chr_len-start)-read_depths[start])/(chr_len - start -1);
		int val = read_depths[start];
		while(val=read_depths[++start]);
//		start++;
//		}
//		else{
	
//		mean= (mean*(chr_len-start)-read_depths[chr_len-1])/(chr_len - 1 - start);
		val = read_depths[chr_len-1];
		while(val == read_depths[--chr_len-1]);
//		}
//		sum = 0;
//		for( i=start;i<chr_len;i++){
//			sum+=(read_depths[i]);
//		}
//		mean = sum / (chr_len-start);
			
		mean= (double)(mean*(chr_len-start)*SLIDING_WINDOW-read_depths[chr_len-1])/(chr_len - 1 - start)/SLIDING_WINDOW;
		sum = 0;
		for( i=start;i<chr_len;i++){
			sum+=(read_depths[i]-mean)*(read_depths[i]-mean);
		}
		std_dev=sqrt((double)sum/(chr_len-start));
	}

	printf("After outlier filtering Chromosome %s mean is %lf\n",snc->chromosome_names[chr_no],mean);

	freeMem(read_depths,sizeof(short)*snc->chromosome_lengths[chr_no]);
	return mean;
}


void calculate_GC_histogram(bam_info *in_bam, sonic *snc, int chr_no){
	int i,j;
	int gc_val;
	int window_per_gc[101];
	long rd_per_gc[101];
	double RD_mean = calculate_chromosome_RD_mean(in_bam,snc,chr_no);
	in_bam->chromosome_mean_rd[chr_no] = RD_mean;
	for( i=0; i<101; i++){
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}
	
	for( i=1; i<snc->chromosome_lengths[chr_no]-SLIDING_WINDOW; i+=SLIDING_WINDOW){
		gc_val = (int) round( sonic_get_gc_content(snc,snc->chromosome_names[chr_no],i,i+SLIDING_WINDOW));

		for(j=0;j<SLIDING_WINDOW;j++){
			rd_per_gc[gc_val] += in_bam->read_depth[chr_no][i+j];
			window_per_gc[gc_val]++;
		}
	}
	in_bam->mean_rd_per_gc[0]= 0.0;
	for( i=1; i<101; i++){
		in_bam->mean_rd_per_gc[i] = READ_LENGTH * (double)rd_per_gc[i] / (window_per_gc[i]);
		if( !isnormal(in_bam->mean_rd_per_gc[i])){
			in_bam->mean_rd_per_gc[i] = 0;
		}
	}
}
*/
double make_global_molecule_mean(short *depths, sonic *snc, int chr){
	
	long bin_count = snc->chromosome_lengths[chr] / MOLECULE_BIN_SIZE;
	double sum = 0;
	int i;
	for(i=0;i<bin_count;i++){
		sum+= depths[i];
	}
	
	return sum/bin_count;
}

double make_global_molecule_std_dev(short *depths, sonic *snc, int chr, double mean){
	double sum = 0;
	long bin_count = snc->chromosome_lengths[chr] / MOLECULE_BIN_SIZE;
	int i;

//	short *depths_copy = malloc(sizeof(short) *bin_count);
//	memcpy(depths_copy,depths,bin_count*sizeof(short));

//	qsort(depths_copy,bin_count,sizeof(short),cmp_short);
	sum = 0;
	for( i=0;i<bin_count;i++){
		sum+=(depths[i]-mean)*(depths[i]-mean);
	}

	double std_dev = sqrt((double)sum/bin_count);

/*
	int start = 0;
	int end = bin_count;
	while( std_dev > ( mean/4)){
	


		int val = depths_copy[start];
		while(val ==depths_copy[++start]);

		val = depths_copy[end-1];
		while(val == depths_copy[--end-1]);

		sum = 0;
		for( i=start;i<end;i++){
			sum+=(depths_copy[i]);
		}
		mean = sum / (end-start);
			
		sum = 0;
		for( i=start;i<end;i++){
			sum+=(depths_copy[i]-mean)*(depths_copy[i]-mean);
		}
		std_dev=sqrt((double)sum/(end-start));
	}
*/
	return std_dev;
}

short *make_molecule_depth_array(vector_t *regions, sonic *snc, int chr){
	long bin_count = 1+snc->chromosome_lengths[chr] / MOLECULE_BIN_SIZE;
	short *depths = malloc(sizeof(short) * bin_count);
	int i, j;

	for( i=0;i<bin_count;i++){
		depths[i] = 0;
	}
	
	for( i=0;i<regions->size;i++){
		interval_10X *molecule = vector_get(regions,i);
		for( j=molecule->start;j<molecule->end;j+=MOLECULE_BIN_SIZE){
			depths[j/MOLECULE_BIN_SIZE]++;	
		}
	}
	return depths;
}

/*
duplication_t *check_cnv(bam_info *in_bam, sonic *snc, int chr, inversion_t *pair){

	double sum;
	double bpAB1_depth;
	double bpAB2_depth;
	double bpCD1_depth;
	double bpCD2_depth;
	int i,j;
	double chr_mean_rd = in_bam->chromosome_mean_rd[chr];

	sum = 0;
	for( i=pair->AB.start1; i<pair->AB.end1;i+=SLIDING_WINDOW){
		double gc_content = in_bam->mean_rd_per_gc[(int)round(sonic_get_gc_content(snc,snc->chromosome_names[chr],i,i+SLIDING_WINDOW))];

		for(j=0;j<SLIDING_WINDOW;j++){
			sum+=in_bam->read_depth[chr][i+j] * (chr_mean_rd/gc_content) ;
		}
	}
	bpAB1_depth = sum/(pair->AB.end1-pair->AB.start1);

	sum = 0;	
	for( i=pair->AB.start2; i<pair->AB.end2;i+=SLIDING_WINDOW){

		double gc_content = in_bam->mean_rd_per_gc[(int)round(sonic_get_gc_content(snc,snc->chromosome_names[chr],i,i+SLIDING_WINDOW))];
		for(j=0;j<SLIDING_WINDOW;j++){
			sum+=in_bam->read_depth[chr][i+j] * (chr_mean_rd/gc_content) ;
		}
	}
	bpAB2_depth = sum/(pair->AB.end2-pair->AB.start2);

	sum = 0;
	for( i=pair->CD.start1; i<pair->CD.end1;i+=SLIDING_WINDOW){
		double gc_content = in_bam->mean_rd_per_gc[(int)round(sonic_get_gc_content(snc,snc->chromosome_names[chr],i,i+SLIDING_WINDOW))];

		for(j=0;j<SLIDING_WINDOW;j++){
			sum+=in_bam->read_depth[chr][i+j] * (chr_mean_rd/gc_content) ;
		}
	}
	bpCD1_depth = sum/(pair->CD.end1-pair->CD.start1);

	sum = 0;	
	for( i=pair->CD.start2; i<pair->CD.end2;i+=SLIDING_WINDOW){
		double gc_content = in_bam->mean_rd_per_gc[(int)round(sonic_get_gc_content(snc,snc->chromosome_names[chr],i,i+SLIDING_WINDOW))];

		for(j=0;j<SLIDING_WINDOW;j++){
			sum+=in_bam->read_depth[chr][i+j] * (chr_mean_rd/gc_content) ;
		}
	}
	bpCD2_depth = sum/(pair->CD.end2-pair->CD.start2);


	duplication_t *new_dup = getMem(sizeof(duplication_t));
	new_dup->break_points = *pair;
	new_dup->bpAB1_depth=bpAB1_depth;
	new_dup->bpAB2_depth=bpAB2_depth;
	new_dup->bpCD1_depth=bpCD1_depth;
	new_dup->bpCD2_depth=bpCD2_depth;
	return new_dup;
}
*/
