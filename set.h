#include "hashtable.h"

#include "common.h"
#include "vector.h"

typedef hashtable_t set_t;

set_t *set_init( size_t table_size, size_t key_size);
vector_t *set_to_vector(set_t *);
//
//  return = { 1    item not in set 
//             0    else            }

int  set_put(set_t *, void *item);
int  set_soft_put(set_t *, void *item);
void set_remove(set_t *, void *item);
int  set_has(set_t *, const void *item);
void set_free(set_t *);
void set_set_remove_function(set_t *, void(*rmv)(void *));
