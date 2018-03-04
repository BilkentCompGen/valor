#include "readbed.h"
#include "valorconfig.h"
#define INITIAL_ARRAY_SIZE 25000



vector_t **read_discordant_bed( char *bed_path){
	int i;
	vector_t **vectors = (vector_t **) getMem(sizeof(vector_t *) *24); //A vector of intervals for each alignment
	fprintf(stderr,"Reading BED file %s.\n", bed_path);

	FILE *bed_file = safe_fopen( bed_path, "r");

	for( i=0; i<24;i++){
		vectors[i] = vector_init(sizeof(interval_discordant),INITIAL_ARRAY_SIZE);
	}
	
	char chr_buf[24];
	int scan_no = 1;
	unsigned long barcode = 0;
	int start1,start2,end1,end2,chr;
	scan_no = fscanf(bed_file,"%s\t%d\t%d\t%d\t%d\t%*c\t%*c",chr_buf,&start1,&end1,&start2,&end2);
	while(  scan_no != -1){

		chr = chr_atoi(chr_buf);
		vector_put(
			vectors[chr],
			&(interval_discordant){start1,end1,start2,end2,barcode}
			);
		scan_no = fscanf(bed_file,"%s\t%d\t%d\t%d\t%d\t%*c\t%*c",chr_buf,&start1,&end1,&start2,&end2);
	}
	fclose(bed_file);
	return vectors;
//	get_library


}


vector_t **read_barcoded_bed( char* bed_path){

	int i;
	vector_t **vectors = (vector_t **) getMem(sizeof(vector_t *) *24); //A vector of intervals for each alignment
	fprintf(stderr,"Reading BED file %s.\n", bed_path);

	FILE *bed_file = safe_fopen( bed_path, "r");

	for( i=0; i<24;i++){
		vectors[i] = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
	}
	
	char chr_buf[24];
	int scan_no = 1;
	unsigned long barcode;
	int start,end,chr;
	scan_no = fscanf(bed_file,"%s\t%d\t%d\t%lu",chr_buf,&start,&end,&barcode);
	while(  scan_no != -1){

		chr = chr_atoi(chr_buf);
		vector_put(
			vectors[chr],
			&(interval_10X){start,end,barcode}
			);
		scan_no = fscanf(bed_file,"%s\t%d\t%d\t%lu",chr_buf,&start,&end,&barcode);
	}
	fclose(bed_file);
	return vectors;
//	get_library_count( in_bam, bam_header->text);

}
vector_t **read_pcs_bed( char* bed_path, int pool_no){

	int i;
	vector_t **vectors = (vector_t **) getMem(sizeof(vector_t *) *24); //A vector of intervals for each alignment
	fprintf(stderr,"Reading BED file %s.\n", bed_path);

	FILE *bed_file = safe_fopen( bed_path, "r");

	for( i=0; i<24;i++){
		vectors[i] = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
	}
	
	char chr_buf[24];
	int scan_no = 1;

	int start,end,chr;
	scan_no = fscanf(bed_file,"%s\t%d\t%d",chr_buf,&start,&end);
	while(  scan_no != -1){

		chr = chr_atoi(chr_buf);
		unsigned long barcode = pool_no;
		vector_put(
			vectors[chr],
			&(interval_10X){start,end,barcode}
			);
		scan_no = fscanf(bed_file,"%s\t%d\t%d",chr_buf,&start,&end);
	}
	fclose(bed_file);
	return vectors;
//	get_library_count( in_bam, bam_header->text);

}

void destroy_beds(vector_t **reads){
	int i;
	for( i = 0; i < 24; i++){
		vector_free(reads[i]);
	}
	freeMem(reads,sizeof(vector_t **));
}
