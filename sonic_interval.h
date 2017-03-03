#ifndef __SONIC_INTERVAL
#define __SONIC_INTERVAL

#include "sonic.h"

typedef enum SONIC_INTERVAL_TYPE {GAP, DUP, REP} sonic_interval_type;

typedef struct _sonic_repeat
{
  char strand;
  char *repeat_type;
  char *repeat_class;
  int repeat_start, repeat_end;
} sonic_repeat;

typedef struct sonic_interval_linked_list
{
  int start;
  int end;
  sonic_repeat *repeat_item;
  struct _interval_linked_list *next;
} sonic_interval_linked_list;

typedef struct sonic_interval_array
{
  int start;
  int end;
  sonic_repeat *repeat_item;
} sonic_interval_array;

typedef struct _sonic_container
{
  int number_of_chromosomes;
  int *number_of_gaps_in_chromosome;
  int *number_of_dups_in_chromosome;
  int *number_of_repeats_in_chromosome;
  sonic_interval_array **gaps;
  sonic_interval_array **dups;
  sonic_interval_array **reps;
} sonic_container;

int sonic_insert(char *, int, int, sonic_interval_type, sonic_repeat *);



#endif
