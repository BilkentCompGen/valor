#ifndef DANG_BIT_SET
#define DANG_BIT_SET
#include "common.h"
#include "vector.h"

typedef struct bit_set{
    size_t size;
    unsigned char *bits;
} bit_set_t;

//static const unsigned char bitmasks[8] = { 1,2,4,8,16,32,64,128};
//static const  unsigned char nbitmasks[8] = {~1,~2,~4,~8,~16,~32,~64,0X7F};
bit_set_t *bit_set_init(size_t);

unsigned char bit_set_get_bit(bit_set_t *, size_t index);

void bit_set_set_bit(bit_set_t *, size_t index, int value);
void bit_set_set_all(bit_set_t *, int value);
void bit_set_flip_bit(bit_set_t *, size_t index);
void bit_set_free(bit_set_t *);
vector_t *bit_set_2_index_vec(bit_set_t *bs);
//TODO
//  bit_set_set_range
//  bit_set_[bitwise_ops]



#endif
