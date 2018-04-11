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

#define PROGRESS_DISTANCE 5001
#define RX_N90_THRESHOLD 4
//What should be the initial size for each array?
#define INITIAL_ARRAY_SIZE 650

//Assumes sorted barcodes
void filter_molecules( vector_t *mols, sonic *snc, int chr){
	unsigned long current_barcode;		
	int i;
	mols->REMOVE_POLICY = REMP_LAZY;

	for(i=0;i<mols->size;i++){
		interval_10X *ival = vector_get(mols,i);
		if(
			sonic_is_gap(snc,snc->chromosome_names[chr]
				,ival->start,ival->end) ||
			sonic_is_satellite(snc,snc->chromosome_names[chr]
				,ival->start,ival->end)){
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
}
vector_t *recover_molecules( vector_t *vector){

	vector_t *regions = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
	vector_t *coverages = vector_init(sizeof(double),INITIAL_ARRAY_SIZE);
	
	printf("Sorting the DNA Intervals\n");
	qsort(vector->items, vector->size, sizeof(void *),barcode_comp);
	int i;

	printf("Done\n");

	int iter = 0;
	double covered;
	printf("Recovering Molecules\n");

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
		for(i=iter+1;i<look_ahead;i++){
			if( I10X_VECTOR_GET(vector,i)->start - merger.end < EXTENSION || merger.start + MOLECULE_EXT >  I10X_VECTOR_GET(vector,i)->start  ){
				merger.end = I10X_VECTOR_GET(vector,i)->end;
				covered+= 180; //interval_size(vector_get(vec tor,i));
			}
			else{
				if( covered / interval_size(&merger) > MIN_COVERAGE &&
						covered / interval_size(&merger) < MAX_COVERAGE
						&& interval_size(&merger) > CLONE_MIN){
					vector_put(regions,&merger);
					vector_put(coverages,&covered);
				}
				merger.start = I10X_VECTOR_GET(vector,i)->start;
				merger.end = I10X_VECTOR_GET(vector,i)->end;
				covered = 180;
			}
		}
		if( covered / interval_size(&merger) > MIN_COVERAGE &&
				covered / interval_size(&merger) < MAX_COVERAGE
				&& interval_size(&merger) > CLONE_MIN){
			vector_put(regions,&merger);
			vector_put(coverages,&covered);
		}
		iter = look_ahead;

		update_progress(iter,vector->size);	

	}
	update_progress(iter,vector->size);		
	printf("\rDone\n");

	//	qsort(regions->items,regions->size,sizeof(void*),interval_start_comp);
#if VALOR_DEBUG
	FILE *write_ptr = fopen(OUT_DIR"/inferredmolecules.out","w+");
	for( i = 0; i < regions->size;i++){
		double *ccp = vector_get(coverages,i);
		fprintf(write_ptr,"%d\t%d\t%d\t%lu\t%lf\n",CUR_CHR,I10X_VECTOR_GET(regions,i)->start,I10X_VECTOR_GET(regions,i)->end,I10X_VECTOR_GET(regions,i)->barcode,*ccp);
	}

	fclose(write_ptr);
#endif
	vector_free(coverages);



	return regions;
}



