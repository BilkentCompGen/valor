#include "vector.h"


void vector_execute_for_all(vector_t *v, void (*foo)(void *)){
        int i;
        for(i=0;i<v->size;i++){
                foo(vector_get(v,i));
        }
}

vector_t *vector_execute_for_all_and_save(vector_t *v, void *(*foo)(void *)){
        int i;
        vector_t *rt = vector_init(VECTOR_VARIABLE_SIZE,v->size);
        for(i=0;i<v->size;i++){
                vector_soft_put(rt,foo(vector_get(v,i)));
        }
        return rt;
}

vector_t *vector_init(size_t item_sizeof, size_t initial_limit){
	vector_t *new_vector = (vector_t *) getMem( sizeof(vector_t));
	new_vector->item_sizeof = item_sizeof;
	new_vector->items = (void **)  getMem( sizeof(void *) * initial_limit);
	new_vector->limit = initial_limit;
	new_vector->size = 0;
	new_vector->REMOVE_POLICY = 0;
	new_vector->fragmental=0;
	new_vector->rmv = &free;
	
	return new_vector;
}

void vector_tabularasa(vector_t *vector){

	int i;
	for(i=0;i<vector->size;i++){
		vector->items[i] = NULL;
	}

	vector->size = 0;
}

void vector_soft_put(vector_t *vector, void *item){
	if(vector->limit == vector->size){
		size_t new_limit = vector->limit +( vector->limit>>1) + 1;
		resizeMem((void **)&(vector->items),vector->limit * sizeof(void *),sizeof(void *) * new_limit);
		vector->limit = new_limit;
	}
	vector->items[vector->size] = item;
	vector->size++;
}

int vector_put(vector_t *vector, void *item){
	if(vector->limit == vector->size){
		size_t new_limit = vector->limit +( vector->limit>>1) + 1;
		resizeMem((void **)&(vector->items),vector->limit * sizeof(void *),sizeof(void *) * new_limit);
		vector->limit = new_limit;
	}
	vector->items[vector->size] = getMem(vector->item_sizeof);

	memcpy( vector->items[vector->size],item,vector->item_sizeof);
	vector->size = vector->size + 1;
	return 0;
}

void vector_update_remove_policy(vector_t *vector, int policy){
	vector->REMOVE_POLICY = policy;
}

int vector_contains(vector_t *vector, void *item){
	int i;
	for(i=0;i<vector->size;i++){
		if(vector->items[i] == NULL){continue;}
		if(memcmp(vector->items[i], item,vector->item_sizeof)==0){
			return i;
		}
	}
	return -1;
}

int vector_comptains(vector_t *vector, void *item, int (*cmp)(const void*, const void*)){
	int i;
	for(i=0;i<vector->size;i++){
		if(cmp(vector->items[i],item)==0){
			return i;
		}
	}
	return -1;
}

int vector_defragment(vector_t *vector){
	if(vector->fragmental==0){
		return 0;
	}
	int i = 0;
	int j = 0;
	while(i<vector->size){
		++i;
		if(vector->items[j]!=NULL){
			++j;
		}
		else{
			vector->fragmental--;
		}
		if(i!=j){
			vector->items[j]=vector->items[i];
	//		vector->items[i]=NULL;
		}
	}
	vector->fragmental=0;
	vector->size = j;
	return 1;
}

void vector_insert(vector_t *vector, void *item, size_t index){
	if(vector->limit == vector->size){
		size_t new_limit = vector->limit +( vector->limit>>1) + 1;
		resizeMem((void **)&(vector->items),vector->limit * sizeof(void *),sizeof(void *) * new_limit);
		vector->limit = new_limit;
	}

	size_t i;
	size_t target = vector->size;
	if(vector->fragmental>0){	
		for(i=index;i<vector->size;i++){
			if(vector->items[i]==NULL){
				target=i;
				vector->fragmental--;
				break;
			}

		}
	}
	for(i=target;i>index;i--){
		vector->items[i]=vector->items[i-1];
	}
	vector->items[index] = getMem(vector->item_sizeof);
	memcpy( vector->items[index],item,vector->item_sizeof);
	vector->size = vector->size + 1;
}

int vector_remove(vector_t *vector, size_t index){
	if(vector->items[index] == NULL || vector->size <= index){
		return -1;
	}

	switch(vector->REMOVE_POLICY){
	case REMP_SORTED:
		vector->rmv(vector->items[index]);

		vector->size--;
		int i;
		for(i=index;i<vector->size;i++){
			vector->items[i] = vector->items[i+1];
		}
	break;
	case REMP_FAST:
		vector->size--;
		vector->rmv(vector->items[index]);
		vector->items[index] = vector->items[vector->size];
	break;
	case REMP_LAZY:
		vector->rmv(vector->items[index]);
		vector->items[index] = NULL;
		vector->fragmental++;
	break;
	default:
		fprintf(stderr,"UNKNOWN POLICY %d\n", vector->REMOVE_POLICY);
		return -2;
	}
	return 0;
}

void vector_clear( vector_t *vector){
	int i;
	for(i=0;i<vector->size;i++){
		vector->rmv(vector->items[i]);
		vector->items[i] = NULL;
	}

	vector->size = 0;
}

//NOT TESTED YET
void vector_zip( vector_t *vector){
	resizeMem((void **)&(vector->items),vector->limit * sizeof(void *),vector->size * sizeof(void *));
}

void *vector_get( vector_t *vector, size_t index){
	#ifdef __DEBUG__
	if(vector->size <= index){
		fprintf(stderr,"Access out of index\n");
		return NULL;
	}
	#endif
	return vector->items[index];
}

void *vector_tail(vector_t *vector){
	return vector->items[vector->size-1];
}

void *vector_head(vector_t *vector){
	return vector->items[0];
}

void vector_free( void *v){
	vector_t *vector = v;
	if(vector==NULL){return;}
	int i;
	for(i = 0; i< vector->size;i++){
		vector->rmv(vector->items[i]);
	}
	freeMem(vector->items,vector->limit * sizeof(void *));
	freeMem(vector,sizeof(vector_t));
}

void vector_swap(vector_t *vec, size_t i1, size_t i2){
	void *tmp = vec->items[i1];
	vec->items[i1] = vec->items[i2];
	vec->items[i2] = tmp;
}

void vector_shuffle(vector_t *vec, unsigned int seed){
	size_t i;
	srand(seed);
	for( i=0;i<vec->size;i++){
		size_t i1 = rand() % vec->size;
		vector_swap(vec,i,i1);
	}
}



int vector_test(int argc, char **argv){
	vector_t *v2 = vector_init(sizeof(int),4);
	int i,j;
	for(j=0;j<10;j++){
//		vector_put(vector,&i);
		vector_put(v2,&j);
	}
	v2->REMOVE_POLICY = REMP_LAZY;

	printf("Vector Size is %zu\n",v2->size);
	vector_remove(v2,3);
	vector_remove(v2,7);
	vector_remove(v2,8);
	i=-26;


	vector_insert(v2,&i,7);
	vector_defragment(v2);

	for(i=0;i<v2->size;i++){
		printf("%d, ",*(int *)vector_get(v2,i));
	}
	printf("\n");

	printf("Vector Size is %zu\n",v2->size);
	for(i=0;i<15;i++){
		printf("Contains: %d -> %d\n",i,vector_contains(v2,&i));
//		printf("%d - %d\n",*((int**)*(void**)v2)[i],v2->items[i]);
	}
	vector_free(v2);
	return 0;
}

void vector_set_remove_function(vector_t *vector, void (*rmv)(void *)){
	vector->rmv = rmv;
}
