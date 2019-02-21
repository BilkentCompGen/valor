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


	sum = 0;
	for( i=0;i<bin_count;i++){
        if(depths[i] < mean / 2.5){
            continue;
        }else if( depths[i] > mean * 2.5){
            continue;
        }
		sum+=(depths[i]-mean)*(depths[i]-mean);
	}

	double std_dev = sqrt((double)sum/bin_count);

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

