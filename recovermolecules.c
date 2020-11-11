#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "progress.h"
#include "recovermolecules.h"
#include "valorconfig.h"
#include "interval10X.h"
#include "vector.h"
#include "common.h"
#include "cmdline.h"

#define PROGRESS_DISTANCE 5001
#define RX_N90_THRESHOLD 2
//What should be the initial size for each array?
#define INITIAL_ARRAY_SIZE 650


vector_t **read_molecules_from_bed(char *filename){
        sonic *snc = sonic_load(NULL);
        vector_t **molecules = malloc(sizeof(vector_t) *snc->number_of_chromosomes);
        int i;
        for(i=0;i<snc->number_of_chromosomes;i++){
                molecules[i] = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
        }

        FILE *bedptr = fopen(filename,"r");
        if(bedptr == NULL){
                fprintf(stderr,"Can't read from %s!\n",filename);
                exit(-1);
        }
        int chr,start,end;
        unsigned long barcode;
        int scan_no;
        scan_no = fscanf(bedptr,"%d\t%d\t%d\t%lu",&chr,&start,&end,&barcode);
        while(scan_no != -1){
                vector_put(molecules[chr],&(interval_10X){start,end,barcode});
                scan_no = fscanf(bedptr,"%d\t%d\t%d\t%lu",&chr,&start,&end,&barcode);
        }
        fclose(bedptr);
        return molecules;
}

void append_molecules_to_bed(vector_t *mols, char *filename, int chr){
	FILE *fptr = fopen(filename,"a+");
	if(fptr==NULL){
		fprintf(stderr,"Can't write to %s!\n",filename);
		exit(-1);
	}
	int i;
	
	for(i=0;i<mols->size;i++){
		interval_10X *val = vector_get(mols,i);
		fprintf(fptr,"%d\t%d\t%d\t%lu\n",chr,val->start,val->end,val->barcode);
	}
	fclose(fptr);
}

//Assumes sorted barcodes
void filter_molecules( vector_t *mols, sonic *snc, int chr){
	unsigned long current_barcode;		
	int i;
	mols->REMOVE_POLICY = REMP_LAZY;
    parameters *params = get_params();
	for(i=0;i<mols->size;i++){
		interval_10X *ival = vector_get(mols,i);
		if(
            (params->filter_gap &&
			sonic_is_gap(snc,snc->chromosome_names[chr]
				,ival->start,ival->end)) ||
            (params->filter_satellite &&
			sonic_is_satellite(snc,snc->chromosome_names[chr]
				,ival->start,ival->end))
            )
        {
			vector_remove(mols,i);
			continue;
		}
	}
	i=0;
	vector_defragment(mols);
	while(i+1<mols->size){
		current_barcode = I10X_VECTOR_GET(mols,i)->barcode;
		if( current_barcode == I10X_VECTOR_GET(mols,i+1)->barcode){
			int j;
			for(j=i+1;j<mols->size
					&& I10X_VECTOR_GET(mols,j)->barcode==current_barcode;j++);
			i=j;

		}
		else{
			vector_remove(mols,i);
			i++;
		}
	}

	vector_defragment(mols);
	mols->REMOVE_POLICY = REMP_SORTED;
}

vector_t *recover_molecules( vector_t *vector){

	vector_t *regions = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
	vector_t *coverages = vector_init(sizeof(double),INITIAL_ARRAY_SIZE);
	
	printf("Sorting the DNA Intervals\n");
	qsort(vector->items, vector->size, sizeof(void *),barcode_comp);
	int i;
	int iter = 0;
	double covered;
	printf("Recovering molecules\n");

	update_progress(iter,vector->size);	

	interval_10X merger = {0,0,0};
	unsigned long current_barcode;

	while( iter < vector->size){
		current_barcode = I10X_VECTOR_GET(vector, iter)->barcode;

		int look_ahead = iter + 1;
		while( look_ahead < vector->size && current_barcode == I10X_VECTOR_GET(vector,look_ahead)->barcode){
			look_ahead++;
		}

		if( look_ahead - iter < RX_N90_THRESHOLD){
			iter = look_ahead;
			continue;
		}
		merger.start = I10X_VECTOR_GET(vector,iter)->start;
		merger.end = I10X_VECTOR_GET(vector,iter)->end;
		merger.barcode = current_barcode;
		covered = 180;
		int read_count = 1;
		for(i=iter+1;i<look_ahead;i++){
			if( I10X_VECTOR_GET(vector,i)->start - merger.end < EXTENSION || merger.start + MOLECULE_EXT >  I10X_VECTOR_GET(vector,i)->start  ){
				merger.end = I10X_VECTOR_GET(vector,i)->end;
				covered+= 180; //interval_size(vector_get(vec tor,i));
				read_count++;
			}
			else{
				if( read_count >= MIN_REQUIRED_READS_IN_MOLECULE &&
						covered / interval_size(&merger) > MIN_COVERAGE &&
						covered / interval_size(&merger) < MAX_COVERAGE
						&& interval_size(&merger) > CLONE_MIN){
					
					vector_put(regions,&merger);
					vector_put(coverages,&covered);
				}
				merger.start = I10X_VECTOR_GET(vector,i)->start;
				merger.end = I10X_VECTOR_GET(vector,i)->end;
				covered = 180;
				read_count = 1;
			}
		}
		if( read_count >= MIN_REQUIRED_READS_IN_MOLECULE &&
			covered / interval_size(&merger) > MIN_COVERAGE &&
				covered / interval_size(&merger) < MAX_COVERAGE
				&& interval_size(&merger) > CLONE_MIN){
			vector_put(regions,&merger);
			vector_put(coverages,&covered);

		}
		iter = look_ahead;

		update_progress(iter,vector->size);	

	}
	update_progress(iter,vector->size);		

	vector_free(coverages);



	return regions;
}



