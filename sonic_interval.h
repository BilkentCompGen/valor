#ifndef __SONIC_INTERVAL
#define __SONIC_INTERVAL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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


int sonic_insert(char *, int, int, sonic_interval_type, struct _sonic_repeat *);

int count_bed_chromosome_entries(sonic_bed_line *, int, char *);
int bed_comp( const void*, const void*);

#endif
