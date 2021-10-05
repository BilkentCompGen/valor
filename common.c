#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "common.h"

// Track memory usage
long long memUsage = 0;


int cmp_int(const void *v1, const void *v2){

    const int64_t *i1 = v1;
    const int64_t *i2 = v2;
    return *i1 - *i2;
}

const char *sv_type_name(sv_type type){
    switch(type){
        case SV_INVERSION:
            return "inversion";
        case SV_DIRECT_DUPLICATION:
            return "duplication";
        case SV_INVERTED_DUPLICATION:
            return "inverted-duplication";
        case SV_TANDEM_DUPLICATION:
            return "tandem-duplication";
        case SV_DELETION:
            return "deletion";
        case SV_TRANSLOCATION:
            return "translocation";
        case SV_INVERTED_TRANSLOCATION:
            return "inverted-translocation";
        case SV_RECIPROCAL:
            return "reciprocal";
        case SV_INVERTED_RECIPROCAL:
            return "inverted-reciprocal";
        default:
            return "unknown";
    }
}



/*
   parameters *init_params(void){

   parameters *params =  get_params();

   params->sonic_file = NULL;
   params->logfile = NULL;
   params->svs_to_find = 0;
   params->low_mem = 0;
   params->chromosome_count = 24;
   params->ploidy = 2;
   (params)->outprefix = NULL;
   (params)->bam_file = NULL;
   (params)->threads = 1;
   params->quasi_clique_lambda = QCLIQUE_LAMBDA;
   params->quasi_clique_gamma = QCLIQUE_GAMMA;
   params->min_qual = MIN_QUAL;
   params->barcode_len = BARCODE_LEN;
   params->sample_size = READ_SAMPLE_SIZE;
   params->max_support = MAX_SUPPORT;
   params->max_frag_size = MAX_FRAG_SIZE;
   params->filter_gap = VALOR_FILTER_GAP;        //true
   params->filter_satellite = VALOR_FILTER_SAT;   //true
   return params;
   }
   */


void print_error( char* msg)
{
    /* print error message and exit */
    fprintf( stderr, "\n%s\n", msg);
    fprintf( stderr, "Invoke parameter -h for help.\n");
    exit( EXIT_COMMON);
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

int identify_read_alignment( bam1_core_t bam_alignment_core, int min, int max)
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

    if(bam_alignment_core.tid != bam_alignment_core.mtid){
        return RPINTER;
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

    if(ret == NULL && new_size != 0){
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

int chr_atoi(char *chromosome){
    if(chromosome[3] > 58){//Some magic number which is between chromosome numbers and letters
        return chromosome[3]-66;
    }
    return atoi(&chromosome[3])-1;
}


int what_is_min_cluster_size(sv_type type,int ploidy){
    double dploidy = ploidy;
    switch(type){
        case SV_DELETION:
            return DELETION_MIN_CLUSTER_SIZE / dploidy;
        case SV_INVERSION:
            return INVERSION_MIN_CLUSTER_SIZE / dploidy;
        case SV_DIRECT_DUPLICATION:
        case SV_INVERTED_DUPLICATION:
            return DUPLICATION_MIN_CLUSTER_SIZE / dploidy;
        case SV_TRANSLOCATION:
        case SV_INVERTED_TRANSLOCATION:
        case SV_RECIPROCAL:
        case SV_INVERTED_RECIPROCAL:
            return TRANSLOCATION_MIN_CLUSTER_SIZE / dploidy;
        case SV_TANDEM_DUPLICATION:
            return TANDEM_DUPLICATION_MIN_CLUSTER_SIZE / dploidy;
        default:
            return -1;
    }
}

// DUP,IDUP,DEL,TRA,INV
sv_type atosv(const char *str){
    if(strcmp(str,"ALL")==0){
        return SV_INVERSION | SV_DIRECT_DUPLICATION | SV_INVERTED_DUPLICATION | SV_DELETION | SV_TRANSLOCATION | SV_INVERTED_TRANSLOCATION;
    }
    if(strcmp(str,"INV")==0){
        return SV_INVERSION;
    }	
    if(strcmp(str,"DUP")==0){
        return SV_DIRECT_DUPLICATION | SV_INVERTED_DUPLICATION;
    }
    if(strcmp(str,"TDUP")==0){
        return SV_TANDEM_DUPLICATION;
    }
    if(strcmp(str,"DEL")==0){
        return SV_DELETION;
    }
    if(strcmp(str,"TRA")==0){
        return SV_TRANSLOCATION | SV_INVERTED_TRANSLOCATION;
    }
    return 0;
}

