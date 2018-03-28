#ifndef __SONIC_INTERVAL
#define __SONIC_INTERVAL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sonic_structures.h"
#include "sonic_reference.h"


int count_bed_chromosome_entries(sonic_bed_line *, int, const char *);
int bed_comp( const void*, const void*);

sonic_interval *sonic_intersect(sonic *, const char *, int, int, sonic_interval_type);
void sonic_print_interval(sonic_interval *);
int sonic_this_interval_intersects(int, int, int, int);
int sonic_is_satellite(sonic *, const char *, int, int);
int sonic_is_gap(sonic *, const char *, int, int);
int sonic_is_segmental_duplication(sonic *, const char *, int, int);
sonic_repeat *sonic_is_mobile_element(sonic *, const char *, int, int, const char *);
float sonic_get_gc_content(sonic *, const char *, int, int);
#endif
