#ifndef __HTABLE
#define __HTABLE
#include <string.h>
#include <stdlib.h>
#include "common.h"

#include "vector.h"

#define INIT_BUCKET_SIZE 4
#define OPTIMAL_LOAD_FACTOR 0.5
#define UPPER_BOUND_LOAD_FACTOR (OPTIMAL_LOAD_FACTOR+0.3)
#define LOWER_BOUND_LOAD_FACTOR (OPTIMAL_LOAD_FACTOR-0.3)

typedef struct __pair_t{
	void *key;
	void *value;
} pair_t;
#define pair_vector_get(V,I) (pair_t *)vector_get((V),(I))

typedef struct __bucket_t{
	void **items;
	size_t item_sizeof;
	size_t limit;
	size_t size;
} bucket_t;


#define BUCKET_VECTOR_TYPE bucket_t

BUCKET_VECTOR_TYPE *bucket_init(size_t item_sizeof, size_t initial_limit);

void bucket_free( void *v);
int bucket_contains(BUCKET_VECTOR_TYPE *vector, void *item);
void *bucket_tail(BUCKET_VECTOR_TYPE *vector);

void bucket_defragment(BUCKET_VECTOR_TYPE *);
int bucket_lazy_remove(BUCKET_VECTOR_TYPE *vector, size_t index);
int bucket_put(BUCKET_VECTOR_TYPE *vector, void *item);
void *bucket_get( BUCKET_VECTOR_TYPE *vector, size_t index);
void bucket_tabularasa(BUCKET_VECTOR_TYPE *vector);
int bucket_remove(BUCKET_VECTOR_TYPE *vector, size_t index);
typedef struct __hashtable_t{
	BUCKET_VECTOR_TYPE **buckets;
	size_t size;
	size_t key_size;
	size_t value_size;
	size_t number_of_items;
	double optimal_load_factor;
	size_t (*hf)(struct __hashtable_t *table, const void *key);
//	void (*key_cpy)(void *, void *, size_t);
//	void (*val_cpy)(void *, void *, size_t);
	int (*key_cmp)(const void *, const void *, size_t);
	void (*key_rmv)(void *);
	void (*val_rmv)(void *);
} hashtable_t;

typedef struct __ht_iter_t{
	hashtable_t *ht;
	size_t bucket_no;
	size_t index;
} ht_iter_t;

ht_iter_t *make_ht_iterator(hashtable_t *table);
int ht_iter_has_next(ht_iter_t *iter);
int ht_iter_next(ht_iter_t *iter);
void *ht_iter_get_key(ht_iter_t *iter);
void *ht_iter_get_value(ht_iter_t *iter);
void ht_iter_free(ht_iter_t *iter);


size_t SuperFastDangHash(hashtable_t *, const void *key);
size_t SuperFastStringHash(hashtable_t *, const void *key);
uint32_t SuperFastHash (const char * data, int len);


vector_t *ht_select_pairs(hashtable_t *,vector_t *);
vector_t *ht_select_pairs_wcmp(hashtable_t *,vector_t *,int (*cmp)(const void *, const void *));
void ht_load_factor_check(hashtable_t *table);
hashtable_t *ht_init( size_t table_size, size_t key_size,size_t item_size);
vector_t *ht_to_vector(hashtable_t *table);
size_t ht_default_hash_function(hashtable_t *table, const void *key);

void *ht_soft_put(hashtable_t *table, void *key);
void *ht_put(hashtable_t *table, const void *key);
void ht_remove(hashtable_t *table, const void *key);

pair_t *ht_get(hashtable_t *table, const void *key);
pair_t *ht_get_wcmp(hashtable_t *table, const void *key, int (*cmp)(const void *, const void *));
void *ht_get_value(hashtable_t *table, const void *key);
int ht_has_key(hashtable_t *table, const void *key);
void ht_free(hashtable_t *table);

void ht_set_key_remove_function(hashtable_t *table, void (*rmv)(void *));
void ht_set_val_remove_function(hashtable_t *table, void (*rmv)(void *));
#endif
