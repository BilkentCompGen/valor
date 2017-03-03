#ifndef __SONIC_INTERVAL
#define __SONIC_INTERVAL

#include "sonic.h"

typedef enum SONIC_ITEM_TYPE {GAP, DUP, REP} sonic_item;

typedef struct _interval_linked_list
{
  int start;
  int end;
  sonic_item item_type;
  void *sonic_repeat;
  struct _interval_linked_list *next;
} interval_linked_list;

typedef struct _sonic_repeat
{
  char strand;
  char *repeat_type;
  char *repeat_class;
  int repeat_start, repeat_end;
} sonic_repeat;



#endif
