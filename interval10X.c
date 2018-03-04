#include "valorconfig.h"
#include "interval10X.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

int interval_overlaps(interval_10X *i1, interval_10X *i2){
	return ((i2->start >= i1->start - MAX_FRAG_SIZE) &&
			(i2->start <= i1->end + MAX_FRAG_SIZE)) ||
		((i2->end >= i1->start - MAX_FRAG_SIZE) &&
		 	(i2->end <= i1->end + MAX_FRAG_SIZE));
}

int interval_pair_overlaps(interval_pair *p1, interval_pair *p2){
	return interval_overlaps( 
	&(interval_10X){p1->start1,p1->end1,p1->barcode},
	&(interval_10X){p2->start1,p2->end1,p2->barcode}
	) && interval_overlaps(
	&(interval_10X){p1->start2,p1->end2,p1->barcode},
	&(interval_10X){p2->start2,p2->end2,p2->barcode}
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
                (i2->start - i1->end >= INV_MIN_SIZE) &&
                (i2->start - i1->end <= INV_MAX_SIZE);
}

