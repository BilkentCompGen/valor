#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"
#include "quotes.h"

// Track memory usage
long long memUsage = 0;

void init_params( parameters** params)
{
	int i;
	/* initialize parameters */
	*params = ( parameters*) getMem( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->sonic = NULL;
	( *params)->reps = NULL;
	( *params)->dups = NULL;
	( *params)->bam_files = NULL;
	( *params)->bam_list_path = NULL;
	( *params)->outprefix = NULL;
	( *params)->bam_file_list = ( char**) getMem( sizeof( char*) * MAX_BAMS);
	( *params)->gaps = NULL;
	( *params)->mei = NULL;
	( *params)->force_read_length = 0;
	( *params)->run_vh = 0; 
	( *params)->run_rd = 0;
	( *params)->run_ns = 0;
	( *params)->run_sr = 0;
	( *params)->threads = 1;
	( *params)->num_bams = 0;
	( *params)->skip_fastq = 0;
	( *params)->skip_sort = 0;
	( *params)->skip_remap = 0;

	for( i = 0; i < MAX_BAMS; i++)
	{
		(*params)->bam_file_list[i] = NULL;
	}
}

void print_params( parameters* params)
{
	int i;
	printf("\n");
	for( i = 0; i < params->num_bams; i++)
	{
		printf( "%-30s%s\n","BAM input:",params->bam_file_list[i]);
		fprintf( logFile,"%-30s%s\n","BAM input:",params->bam_file_list[i]);
	}
	printf( "%-30s%s\n","Reference genome:", params->ref_genome);
	printf( "%-30s%s\n","Repeat coordinates:", params->reps);
	printf( "%-30s%s\n","Duplication coorinates:", params->dups);
	printf( "%-30s%s\n","Assembly Gap coordinates:", params->gaps);
	printf( "%-30s%s\n","Mobile Elements:", params->mei);

	fprintf( logFile, "%-30s%s\n","Reference genome:", params->ref_genome);
	fprintf( logFile, "%-30s%s\n","Repeat coordinates:", params->reps);
	fprintf( logFile, "%-30s%s\n","Duplication coorinates:", params->dups);
	fprintf( logFile, "%-30s%s\n","Assembly Gap coordinates:", params->gaps);
	fprintf( logFile, "%-30s%s\n","Mobile Elements:", params->mei);
}

void print_error( char* msg)
{
	/* print error message and exit */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	exit( EXIT_COMMON);
}

void print_quote( void)
{
	/* print a quote by the Doctor */

	int quotenum;

	srand(time(NULL));
	quotenum = rand() % NUM_QUOTES;
	fprintf( stderr, "\n\t%s\n\n", quotes[quotenum]);
}

FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[VALOR INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);

	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	gzFile file;
	char err[500];

	file = gzopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[VALOR INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);		
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];

	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[VALOR INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}


int is_proper( int flag)
{
	if ( (flag & BAM_FPAIRED) != 0 && (flag & BAM_FSECONDARY) == 0 && (flag & BAM_FSUPPLEMENTARY) == 0 && (flag & BAM_FDUP) == 0 && (flag & BAM_FQCFAIL) == 0)
		return 1;

	return 0;
}



int is_discordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPAIRED) == 0) 
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPSEND;
	}

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return RPINV;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return RPINV;
	}

	if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
	{
		/* Read is placed BEFORE its mate */
		if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
		{
			/* -+ orientation = tandem duplication */
			return RPTDUP;
		}
	}
	else
	{
		/* Read is placed AFTER its mate */
		if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
		{
			/* +- orientation = tandem duplication */
			return RPTDUP;
		}
	}

	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if(abs(bam_alignment_core.isize) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}


int is_alt_concordant( int p1, int p2, int flag, char strand1, char strand2, int min, int max){
	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPSEND;
	}
	/*
	if( ( flag & BAM_FPROPER_PAIR) == 0)
	{
		//Not proper pair
		return RPUNMAPPED;
	}*/

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( strand1 != 0 && strand2 != 0)
	{
		/* -- orientation = inversion */
		return RPMM;
	}

	if( strand1 == 0 && strand2 == 0)
	{
		/* ++ orientation = inversion */
		return RPPP;
	}

	if( abs(p1-p2) > DUP_MIN_DIST){
//		if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
		{
			/* Read is placed BEFORE its mate */
			if(  strand1 != 0 && strand2 == 0)
			{
				/* -+ orientation = tandem duplication */
				return RPTDUPPM; //was 0 before
			}
		}
//		else
		{
			/* Read is placed AFTER its mate */
			if(  strand1 == 0 && strand2 != 0)
			{
				/* +- orientation = tandem duplication */
				return RPTDUPMP; //was 0 before
			}
		}


	}
	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(p1-p2) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if(abs(p1-p2) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}



int is_concordant( bam1_core_t bam_alignment_core, int min, int max)
{
	int flag = bam_alignment_core.flag;

	if( ( flag & BAM_FPAIRED) == 0)
	{
		/* Read is single-end. Skip this by calling it concordant */
		return RPSEND;
	}
	/*
	if( ( flag & BAM_FPROPER_PAIR) == 0)
	{
		//Not proper pair
		return RPUNMAPPED;
	}*/

	if( ( flag & BAM_FUNMAP) != 0)  // c.a.
	{
		/* Read unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FMUNMAP) != 0) // c.a.
	{
		/* Mate unmapped; Orphan or OEA */
		return RPUNMAPPED;
	}

	if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) != 0)
	{
		/* -- orientation = inversion */
		return RPMM;
	}

	if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) == 0)
	{
		/* ++ orientation = inversion */
		return RPPP;
	}

	if( bam_alignment_core.tid != bam_alignment_core.mtid)
	{
		/* On different chromosomes */
		return RPUNMAPPED;
	}

	if( abs(bam_alignment_core.pos-bam_alignment_core.mpos) > DUP_MIN_DIST){
//		if( bam_alignment_core.pos <= bam_alignment_core.mpos) // c.a.
		{
			/* Read is placed BEFORE its mate */
			if( ( flag & BAM_FREVERSE) != 0 && ( flag & BAM_FMREVERSE) == 0)
			{
				/* -+ orientation = tandem duplication */
				return RPTDUPPM; //was 0 before
			}
		}
//		else
		{
			/* Read is placed AFTER its mate */
			if( ( flag & BAM_FREVERSE) == 0 && ( flag & BAM_FMREVERSE) != 0)
			{
				/* +- orientation = tandem duplication */
				return RPTDUPMP; //was 0 before
			}
		}


	}
	/* Passed all of the above. proper pair, both mapped, in +- orientation. Now check the isize */
	if( abs(bam_alignment_core.isize) < min) // c.a.
	{
		/* Deletion or Insertion */
		return RPINS;
	}
	else if(abs(bam_alignment_core.isize) > max)
	{
		return RPDEL;
	}

	/* All passed. Read is concordant */
	return RPCONC;
}

/* Decode 4-bit encoded bases to their corresponding characters */
char base_as_char( int base_as_int)
{
	if( base_as_int == 1)
	{
		return 'A';
	}
	else if( base_as_int == 2)
	{
		return 'C';
	}
	else if( base_as_int == 4)
	{
		return 'G';
	}
	else if( base_as_int == 8)
	{
		return 'T';
	}
	else if( base_as_int == 15)
	{
		return 'N';
	}
	return 0;
}

/* Return the complement of a base */
char complement_char( char base)
{
	switch( base)
	{
	case 'A':
		return 'T';
		break;
	case 'C':
		return 'G';
		break;
	case 'G':
		return 'C';
		break;
	case 'T':
		return 'A';
		break;
	default:
		return 'N';
		break;
	}
	return 'X';
}

/* Add 33 to the interger value of the qual characters to convert them to ASCII */
void qual_to_ascii( char* qual)
{
	int i;
	for( i = 0; i < strlen( qual); i++)
	{
		qual[i] = qual[i] + 33;
	}
}

/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	if (source != NULL)
	{
		( *target) = ( char*) getMem( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
		( *target) = NULL;
	}
}


/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}

int compare_size_int( const void* p, const void* q)
{
	int i = *( const int*) p;
	int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void* getMem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( size / 1048576.0));
		exit( 0);
	}

	memUsage = memUsage + size;
	return ret;
}

void resizeMem( void **ptr, size_t old_size, size_t new_size){
	void *ret;

	ret = realloc(*ptr,new_size);

	if(ret == NULL){
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( (new_size-old_size) / 1048576.0));
		exit( 0);

	}
	*ptr = ret;	
	memUsage = memUsage + new_size - old_size;
}

void freeMem( void* ptr, size_t size)
{
	memUsage = memUsage - size;
	free( ptr);
}

double getMemUsage()
{
	return memUsage / 1048576.0;
}

/*
int count_bed_lines(FILE *bed_file)
{
	int number_of_lines;
	char line[MAX_LENGTH];
	char *return_value;

	number_of_lines = 0;
	while (!feof(bed_file)){
		return_value = fgets(line, MAX_LENGTH, bed_file);
		if (feof(bed_file))
			break;
		if (line[0] != 0)
			number_of_lines++;
	}

	return number_of_lines;
}
*/
int findChroIndex(ref_genome* ref, char* chroName)
{
	int i;
	for(i=0;i<ref->chrom_count;i++)
	{
		if(!strcmp(chroName, ref->chrom_names[i]))
			return i;
	}
	return -1;
}


int chr_atoi(char *chromosome){
        if(chromosome[3] > 58){//Some magic number which is between chromosome numbers and letters
                return chromosome[3]-66;
        }
        return atoi(&chromosome[3])-1;
}
