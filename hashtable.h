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

typedef struct __hashtable_t{
	vector_t **buckets;
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
void *ht_put(hashtable_t *table, void *key);
void ht_remove(hashtable_t *table, void *key);

pair_t *ht_get(hashtable_t *table, void *key);
pair_t *ht_get_wcmp(hashtable_t *table, void *key, int (*cmp)(const void *, const void *));
void *ht_get_value(hashtable_t *table, void *key);
int ht_has_key(hashtable_t *table, void *key);
void ht_free(hashtable_t *table);

void ht_set_key_remove_function(hashtable_t *table, void (*rmv)(void *));
void ht_set_val_remove_function(hashtable_t *table, void (*rmv)(void *));
#endif
