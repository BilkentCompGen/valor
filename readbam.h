#ifndef __READBAM
#define __READBAM

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "interval10X.h"
#include "vector.h"
#include "sonic.h"
typedef struct __bam_vector_pack{
	vector_t *concordants;
	vector_t *mm_discordants;
	vector_t *pp_discordants;
	vector_t *pm_discordants;
	vector_t *mp_discordants;
}bam_vector_pack;


typedef struct simple_interval{
	int tid;
	int start;
	int end;
	char strand;
}simple_interval;


typedef struct alt_read{
	char *read_name;
	unsigned long barcode;
	simple_interval *positions;
	short flag;
	char count;
}alt_read;

void free_alt_read(void *read);
int altcomp(const void *, const void *);
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
} bam_info;
vector_t *dang_string_tokenize(const char *str, const  char *delimiter);
bam_stats *calculate_bam_statistics(bam_info *, char *bam_path, int number_of_reads);
bam_vector_pack *make_bam_vector_pack();
int read_pcs_bam( bam_info* in_bam, char *bam_path, int pool_no, bam_vector_pack *pack);
bam_vector_pack **read_10X_bam( bam_info* in_bam, char *bam_path, sonic *snc);
void destroy_bams( bam_vector_pack* reads);
#endif

