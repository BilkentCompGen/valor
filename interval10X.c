#include "valorconfig.h"
#include "interval10X.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cgranges/cgranges.h"



int64_t *union_overlaps( int64_t *a, int64_t lena, int64_t *b, int64_t lenb, int64_t *lenret){

    int64_t ai = 0;
    int64_t bi = 0;
    int64_t index = 0;

    int64_t *ret = malloc((lena + lenb) * sizeof(int64_t));
    
    while( ai < lena && bi < lenb){
        if(a[ai] == b[bi]){
            ret[index] = a[ai];
            ++index;
            ++ai;
            ++bi;
        }
        else if( a[ai] > b[bi]){
            ret[index] = b[bi];
            ++index;
            ++bi;
        }
        else{
            ret[index] = a[ai];
            ++index;
            ++ai;
        }

    }
    *lenret = index;
    return ret;
}
int64_t *intersection_overlaps( int64_t *a, int64_t lena, int64_t *b, int64_t lenb, int64_t *lenret){

    int64_t ai = 0;
    int64_t bi = 0;
    int64_t index = 0;

    int64_t *ret = malloc(MAX(lena, lenb) * sizeof(int64_t));
    
    while( ai < lena && bi < lenb){
        if(a[ai] == b[bi]){
            ret[index] = a[ai];
            ++index;
            ++ai;
            ++bi;
        }
        else if( a[ai] > b[bi]){
            ++bi;
        }
        else{
            ++ai;
        }

    }
    *lenret = index;
    return ret;
}

void replace_cg_labels( cgranges_t *cg, int64_t *arr, int64_t len){
    int i;

    for(i=0;i< len;++i){
        arr[i] = cr_label(cg,arr[i]);
    }

}


int interval_start_comp(const void *v1, const void *v2){
        interval_10X *i1 = *(interval_10X **)v1;
        interval_10X *i2 = *(interval_10X **)v2;
	return i1->start-i2->start;
}
 
int interval_comp(const void *v1,const void *v2){
        interval_10X *i1 = *(interval_10X **)v1;
        interval_10X *i2 = *(interval_10X **)v2;
        if (i1->start < i2->start) return -1;
        if (i2->start < i1->start) return 1;
        if (i1->end < i2->end) return -1;
        if (i2->end < i1->end) return 1;
        return 0;
}

int in_range(int a, int b, int range){
    return (a-range<b) && (a+range >b);
}

//Probably this will achieve a stability by using interval_comp in case of equality
int barcode_comp(const void *v1, const void *v2){
	interval_10X *i1 = (*(void **)v1);
	interval_10X *i2 = (*(void **)v2);
	if ( i1->barcode == i2->barcode){
		return interval_comp(v1,v2);
	}
	return i1->barcode > i2->barcode?1:-1;
}

int interval_pair_comp(const void *v1, const void *v2){
	interval_pair *p1 = *(void **)v1;
	interval_pair *p2 = *(void **)v2;
	if(p1->start1 - p2->start1 == 0){
		return p1->start2 - p2->start2;
	}
	return p1->start1 - p2->start1;
}

int discordant_barcode_comp(const void *v1, const void *v2){
	interval_pair *i1 = (*(void **)v1);
	interval_pair *i2 = (*(void **)v2);
	if ( i1->barcode == i2->barcode){
		return interval_pair_comp(v1,v2);
	}
	return i1->barcode > i2->barcode?1:-1;
}
//Distance for inversions
int i_distance(int start1, int start2, int end1, int end2){
	if(start1 < start2){
		return i_distance(start2,start1,end2,end1);
	}
	return start1-end2;
}

//Distance for duplications
int d_distance(int start1, int start2, int end1, int end2){
	if(start1 < start2){
		return d_distance(start2, start1, end2, end1);
	}
	return end1-start2;
}

int interval_outer_distance(interval_10X a, interval_10X b){
	if( a.end > b.end){
		return a.end - b.start;
	}
	return b.end - a.start;
}
int interval_inner_distance(interval_10X a, interval_10X b){
	if( a.start < b.start){
		return b.start - a.end;
	}
	return a.start - b.end;
}

//Fix this for varying length
char *decode_ten_x_barcode( unsigned long barcode){
	int i;
	char *text = malloc(17); 
	unsigned long mask = 3;
	char *ACGT = "ACGT";
	for(i=0;i<16;i++){
		unsigned long res = barcode & mask;
		mask= mask << 2;
		text[i]=ACGT[res];
	}	
	text[16]=0;
	return text;
}

/*
 *	Borrowed from tardis
 *	char* to unsigned long barcode encoder
 */
unsigned long encode_ten_x_barcode(unsigned char* source){
	int i, len;
	unsigned long next_digit, result;

	if (source == NULL){
		return -1;
	}
	result = 0;
	next_digit = 0;
	len = strlen((char *)source);
	for (i = 0; i < len; i++){
		switch(source[i]){
		case 'A':
			next_digit = 0;
			result = (result << 2)|next_digit;
			break;
		case 'C':
			next_digit = 1;
			result = (result << 2)|next_digit;
			break;
		case 'G' :
			next_digit = 2;
			result = (result << 2)|next_digit;
			break;
		case 'T':
			next_digit = 3;
			result = (result << 2)|next_digit;
			break;
		default:
			break;
		}
	}
	return result;
}

int interval_overlaps(interval_10X *i1, interval_10X *i2, int relax){
	return ((i2->start >= i1->start - relax) &&
			(i2->start <= i1->end + relax)) ||
		((i2->end >= i1->start - relax) &&
		 	(i2->end <= i1->end + relax));
}

int interval_pair_overlaps(interval_pair *p1, interval_pair *p2,int relax){
	return interval_overlaps( 
	&(interval_10X){p1->start1,p1->end1,p1->barcode},
	&(interval_10X){p2->start1,p2->end1,p2->barcode},
	relax
	) && interval_overlaps(
	&(interval_10X){p1->start2,p1->end2,p1->barcode},
	&(interval_10X){p2->start2,p2->end2,p2->barcode},
	relax
	);
}

void interval_pair_intersect( interval_pair *p1, interval_pair *p2){
	p1->start1 = MAX(p1->start1,p2->start1);
	p1->start2 = MAX(p1->start2,p2->start2);
	p1->end1 = MIN(p1->end1,p2->end1);
	p1->end2 = MIN(p1->end2,p2->end2);
}

int interval_size(interval_10X *this){
	return this->end - this->start +1;
}

int interval_distance(interval_10X *this, interval_10X *that){
	if(this->start < that->start){
		return interval_distance(that,this);
	}
	return this->start - that->end;
}

//Assumes intervals are from same chromosome
int interval_can_pair(interval_10X *i1, interval_10X *i2){
	        return (i1->barcode == i2->barcode) &&
                (interval_size(i1) + interval_size(i2) - 2 >= CLONE_MIN) &&
                (interval_size(i1) + interval_size(i2) - 2 <= CLONE_MAX) &&
                (i2->start - i1->end >= CLONE_MIN_DIST) &&
                (i2->start - i1->end <= CLONE_MAX_DIST);
}



#define SCL_INIT_LIMIT 10000


vector_t *discover_split_molecules(vector_t *regions){
	int i,j;
	vector_t *smolecules = vector_init(sizeof(splitmolecule_t),SCL_INIT_LIMIT);
	interval_10X *ii;
	interval_10X *ij;

	for(i=0;i<regions->size;i++){
		for(j=i+1;j<regions->size;j++){
			ii = vector_get(regions,i);
			ij = vector_get(regions,j);
			if( interval_can_pair(ii,ij)){
				splitmolecule_t *new_molecule = splitmolecule_init(ii,ij);
				vector_soft_put(smolecules,new_molecule);
			}else if( interval_can_pair(ij,ii)){
				splitmolecule_t *new_molecule = splitmolecule_init(ij,ii);
				vector_soft_put(smolecules,new_molecule);
			}
			//	VALOR_LOG("%d-%d\t%d-%d\n",i,j,ii->barcode,ij->barcode);
			if( ij->barcode != ii->barcode){
				break;
			}	
		}

	}

	return smolecules;
}


splitmolecule_t *splitmolecule_copy(splitmolecule_t *scl){
	splitmolecule_t *ncl = getMem(sizeof(splitmolecule_t));

	ncl->start1 = scl->start1;
	ncl->start2 = scl->start2;
	ncl->end1 = scl->end1;
	ncl->end2 = scl->end2;
	ncl->barcode = scl->barcode;
	return ncl;
}
size_t scl_binary_search(vector_t *intervals, splitmolecule_t *key){
	if(intervals->size == 0){return -1;}
	long first, last;
	long mid = 0;
	first =0;
	last = intervals->size - 1;
	int counter = 0;	
	while( first < last){
		mid = (first + last)/2;
		if(IDIS_VECTOR_GET(intervals,mid)->end1 < key->start1){
			first = mid + 1;
		}
		else{
			last = mid - 1;
		}
		counter ++;
	}
	return mid;
}
void splitmolecule_destroy(splitmolecule_t *molecule){
	free(molecule);
}
splitmolecule_t *splitmolecule_init(interval_10X *i1, interval_10X *i2){
	if( i1->barcode != i2->barcode){
		fprintf(stderr,"Can't create split molecule of different barcodes\n");
		return NULL;
	}
	splitmolecule_t *new_molecule = malloc(sizeof(splitmolecule_t));
	new_molecule->start1 = i1->start;
	new_molecule->start2 = i2->start;
	new_molecule->end1 = i1->end;
	new_molecule->end2 = i2->end;
	new_molecule->barcode = i1->barcode;
	return new_molecule;
}
