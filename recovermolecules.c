/*
	TODO:			
	Learn to use is_concordant from common.h
*/
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

//int interval_size(interval *this){
//        return this->end - this->start +1;
//}
#define PROGRESS_DISTANCE 5001
#define RX_N90_THRESHOLD 4
//What should be the initial size for each array?
#define INITIAL_ARRAY_SIZE 650

//Assumes sorted barcodes
void filter_molecules( vector_t *mols){
	unsigned long current_barcode;		
	int i = 0;
	mols->REMOVE_POLICY = REMP_LAZY;

	while(i+1<mols->size){
		current_barcode = I10X_VECTOR_GET(mols,i)->barcode;
		if( current_barcode == I10X_VECTOR_GET(mols,i+1)->barcode){
			int j;
			for(j=i+1;j<mols->size
				 && I10X_VECTOR_GET(mols,j)->barcode==current_barcode;j++);
			i=j;

		}
		else{
//			VALOR_LOG("%lu - %lu are not equal at %d\n",current_barcode,I10X_VECTOR_GET(mols,i+1)->barcode,i);	
	//		printf("removed %d: current barcode %lu: next barcode %lu\n",i,current_barcode,I10X_VECTOR_GET(mols,1+i)->barcode);
			vector_remove(mols,i);
			i++;
		}
	}
	int j=0;
	for(i=0;i<mols->size;i++){
		if(vector_get(mols,i)!=NULL){
			j++;
		}
	}
//	printf("vsize %zu\n",j);
	vector_defragment(mols);
//	printf("vsize %zu\n",mols->size);

	//qsort(mols->items,mols->size,sizeof(void*),interval_start_comp);
}
vector_t *recover_molecules( vector_t *vector){
	printf("Sorting the DNA Intervals\n");
	vector_t *regions = vector_init(sizeof(interval_10X),INITIAL_ARRAY_SIZE);
	vector_t *coverages = vector_init(sizeof(double),INITIAL_ARRAY_SIZE);
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
	#if DEVELOPMENT_
	FILE *write_ptr = fopen(OUT_DIR"inferredmolecules.out","a+");
	for( i = 0; i < regions->size;i++){
		double *ccp = vector_get(coverages,i);
		fprintf(write_ptr,"%d\t%d\t%d\t%lu\t%lf\n",CUR_CHR,I10X_VECTOR_GET(regions,i)->start,I10X_VECTOR_GET(regions,i)->end,I10X_VECTOR_GET(regions,i)->barcode,*ccp);
	}

	fclose(write_ptr);
	#endif
	vector_free(coverages);
	


	return regions;
}



