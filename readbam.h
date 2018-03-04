#ifndef __READBAM
#define __READBAM

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "interval10X.h"
#include "vector.h"
#include "sonic.h"
typedef struct __bam_vector_pack{
	vector_t **concordants;
	vector_t **mm_discordants;
	vector_t **pp_discordants;
	vector_t **pm_dup;
	vector_t **mp_dup;
	int count;
}bam_vector_pack;

typedef struct __bam_stats{
	double read_length_mean;
	double read_length_std_dev;
} bam_stats;

enum gender {male,female};

typedef struct _bam_info
{
        htsFile* bam_file; /* file pointer to the BAM file */
        enum gender sample_gender; /* gender of the sample */
        char* sample_name; /* name of the sample, parsed from SM in the BAM header */
        int num_libraries; /* number of libraries, counted from the RG tags in the BAM header */
        struct library_properties** libraries; /* each library_properties struct holds statistical/other info */
        int read_count;
        float mean_rd_per_gc[101]; // Vector of GC percentage arrays of each chromosome
        double *chromosome_mean_rd;
        short **read_depth; // Vector of short arrays for read_depth of each window
        int number_of_chr;
        char **chromosome_names;
} bam_info;
bam_stats *calculate_bam_statistics(bam_info *, char *bam_path, int number_of_reads);
bam_vector_pack *make_bam_vector_pack();
int read_pcs_bam( bam_info* in_bam, char *bam_path, int pool_no, bam_vector_pack *pack);
bam_vector_pack *read_10X_bam( bam_info* in_bam, char *bam_path, sonic *snc);
void destroy_bams( bam_vector_pack* reads);
#endif

