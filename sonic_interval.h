#ifndef __SONIC_INTERVAL
#define __SONIC_INTERVAL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sonic_structures.h"
#include "sonic_reference.h"


int count_bed_chromosome_entries(sonic_bed_line *, int, char *);
int bed_comp( const void*, const void*);

sonic_interval *sonic_intersect(sonic *, char *, int, int, sonic_interval_type);
void sonic_print_interval(sonic_interval *);
int sonic_this_interval_intersects(int, int, int, int);
int sonic_is_satellite(sonic *, char *, int, int);
int sonic_is_gap(sonic *, char *, int, int);
int sonic_is_segmental_duplication(sonic *, char *, int, int);
sonic_repeat *sonic_is_mobile_element(sonic *, char *, int, int, char *);

#endif
