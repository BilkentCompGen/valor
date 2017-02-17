#ifndef __SONIC
#define __SONIC

#define SONIC_MAGIC 42

#include <zlib.h>
#include "common.h"
#include <string.h>
#include <stdlib.h>


int make_sonic(parameters *);
int load_sonic(parameters *);


#endif
