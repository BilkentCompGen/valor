#include "statistics.h"
double molecule_mean(vector_t *regions){
        int i;
        double sum = 0;
        for(i=0;i<regions->size;i++){
                interval_10X *tmp = vector_get(regions,i);
                sum+=interval_size(tmp);
        }
        return sum/regions->size;
}

double molecule_std(vector_t *regions, double mean){
        int i;
        double var = 0;
        for(i=0;i<regions->size;i++){
                interval_10X *tmp = vector_get(regions,i);
                int isize = interval_size(tmp);
                var+=(isize-mean)*(isize-mean);
        }
        return sqrt(var/regions->size);
}

vector_t *group_overlapping_molecules(vector_t *molecules){
	if(molecules->size==0){ return vector_init(sizeof(vector_t),1);}
	int i;
	int index = 0;
	vector_t *vector_of_groups = vector_init(sizeof(vector_t),16);
	vector_of_groups->rmv = &vector_free;
	vector_soft_put(vector_of_groups,
			vector_init(sizeof(interval_10X),16));	
	interval_10X *current = vector_get(molecules,0);

	vector_put(vector_get(vector_of_groups,index),current);

	for(i=1;i<molecules->size;i++){
		if(!interval_overlaps(current,vector_get(molecules,i))){
			index++;
			vector_soft_put(vector_of_groups,
				vector_init(sizeof(interval_10X),16));
			current=vector_get(molecules,i);
			vector_put(vector_get(vector_of_groups,index),
				current);
		}
		else{
			vector_put(vector_get(vector_of_groups,index),
				vector_get(molecules,i));
		}
	}
	return vector_of_groups;
}

double molecule_group_std(vector_t *groups, double mean){
	int i;
	double group_means[groups->size];
	double var = 0;
	for(i=0;i<groups->size;i++){
		group_means[i]=molecule_mean(vector_get(groups,i));
	}

	for(i=0;i<groups->size;i++){
		double isize = group_means[i];
                var+=(isize-mean)*(isize-mean);
        }
        return sqrt(var/groups->size);

}
