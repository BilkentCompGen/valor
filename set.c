#include "set.h"

set_t *set_init( size_t table_size, size_t key_size){
    return  ht_init(table_size, key_size, 0);
}

vector_t *set_to_vector(set_t *set){
    int i,j;
    vector_t *v = vector_init(sizeof(set->key_size),set->number_of_items);
    for(i=0;i<set->size;i++){
        bucket_t *bucket = set->buckets[i];
        for(j=0;j<bucket->size;j++){
            pair_t *dummy = bucket_get(bucket,j);
            vector_put(v,dummy->key);
        }
    }
    return v;
}

int set_put(set_t *set, void *item){
    if( set_has(set, item)){
        return 0;
    }
    ht_put(set,item);
    return 1;
}

void set_remove(set_t *set, void *item){
    ht_remove(set,item);
}

int set_has(set_t *set, void *item){
    return ht_has_key(set,item);
}

void set_free(set_t *set){
    ht_free(set);
}

void set_set_remove_function(set_t *set, void (*rmv)(void *)){
    set->key_rmv = rmv;
}
