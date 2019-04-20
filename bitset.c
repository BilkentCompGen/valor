#include "bitset.h"


static const unsigned char bitmasks[8] = { 1,2,4,8,16,32,64,128};
static const unsigned char nbitmasks[8] = {~1,~2,~4,~8,~16,~32,~64,0X7F};

bit_set_t *bit_set_init(size_t bit_count){
    bit_set_t *new_bits = getMem(sizeof(bit_set_t));
    new_bits->size = bit_count;
    size_t byte_count = bit_count / 8 + (bit_count % 8 !=0?1:0);
    int more = sizeof(long) - ( byte_count % sizeof(long));
    byte_count+=more;
    new_bits->bits = calloc(sizeof(unsigned char), byte_count);
    return new_bits;
}

bit_set_t __bit_set_make(size_t bit_count){
    bit_set_t new_bits;
    new_bits.size = bit_count;
    size_t byte_count = bit_count / 8 + (bit_count % 8 !=0?1:0);
    int more = sizeof(long) - ( byte_count % sizeof(long));
    byte_count+=more;
    new_bits.bits = calloc(sizeof(unsigned char), byte_count);
    return new_bits;
}

unsigned char bit_set_get_bit(bit_set_t *bs, size_t bitdex){
    size_t index = bitdex / 8;
    int bit_index = bitdex % 8;
    if(index >= bs->size){return 0;}	    
    return (bs->bits[index] & bitmasks[bit_index]) >> (bit_index) ;
}

void bit_set_set_bit( bit_set_t *bs, size_t bitdex, int val){
    size_t index = bitdex / 8;
    unsigned char bit_index = bitdex % 8;
    if(val){
        bs->bits[index] = bs->bits[index] | bitmasks[bit_index];
    }
    else{
        bs->bits[index] = bs->bits[index] & nbitmasks[bit_index];
    }
}


void bit_set_set_all( bit_set_t *bs, int value){
    int i;
    size_t llimit = bs->size / (8*sizeof(long)) + 1;

    unsigned long *longs = (unsigned long *) bs->bits;
    for(i=0;i<llimit;i++){
        longs[i] = (value?-1:0);
    }
}

void bit_set_flip_bit(bit_set_t *bs, size_t bitdex){
    size_t index = bitdex / 8;
    unsigned char bit_index = bitdex % 8;
    int val = bs->bits[index] & bitmasks[bit_index];
    if(!val){
        bs->bits[index] = bs->bits[index] | bitmasks[bit_index];
    }
    else{
        bs->bits[index] = bs->bits[index] & nbitmasks[bit_index];
    }
}

vector_t *bit_set_2_index_vec(bit_set_t *bs){
	int i;
	vector_t *out = vector_init(sizeof(int),bs->size);
	
	for(i=0;i<bs->size;i++){
		if(bit_set_get_bit(bs,i)){
			vector_put(out,&i);
		}	
	}
	return out;
}

void bit_set_free(bit_set_t *bs){
    free(bs->bits);
    free(bs);
}
