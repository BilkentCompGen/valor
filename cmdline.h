#ifndef __COMMANDLINE
#define __COMMANDLINE

#include "sonic/sonic.h"
#include "common.h"

int parse_command_line( int, char**, parameters*, sonic *);
void parse_bam_list( parameters** params);
void print_help( void);

#endif
