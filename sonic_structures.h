#ifndef __SONIC_STRUCTURES
#define __SONIC_STRUCTURES


typedef enum SONIC_INTERVAL_TYPE {SONIC_GAP, SONIC_DUP, SONIC_REP} sonic_interval_type;

typedef struct _sonic_repeat
{
  char strand;
  char *repeat_type;
  char *repeat_class;
  int repeat_start, repeat_end;
  int mei_code;
} sonic_repeat;


typedef struct sonic_interval_linked_list
{
  int start;
  int end;
  struct _sonic_repeat *repeat_item;
  struct _interval_linked_list *next;
} sonic_interval_linked_list;

typedef struct _sonic_interval
{
  int start;
  int end;
  struct _sonic_repeat *repeat_item;
} sonic_interval;


typedef struct _sonic_bed_line
{
  char chromosome[255];
  int start;
  int end;
  struct _sonic_repeat *repeat_item;
} sonic_bed_line;

typedef struct _sonic
{
  int number_of_chromosomes;
  int *number_of_gaps_in_chromosome;
  int *number_of_dups_in_chromosome;
  int *number_of_repeats_in_chromosome;
  int *chromosome_lengths;
  long genome_length;
  char **chromosome_gc_profile;
  char **chromosome_names;
  int last_chromosome_index;
  struct _sonic_interval **gaps;
  struct _sonic_interval **dups;
  struct _sonic_interval **reps;
} sonic;


#endif
