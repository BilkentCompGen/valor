#ifndef __READBAM
#define __READBAM

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "bitset.h"
#include "interval10X.h"
#include "vector.h"
#include "sonic.h"
typedef struct __bam_vector_pack{
	vector_t *concordants;
	vector_t *mm_discordants;
	vector_t *pp_discordants;
	vector_t *pm_discordants;
	vector_t *mp_discordants;

	vector_t *inter_pm;
	vector_t *inter_mp;
	vector_t *inter_pp;
	vector_t *inter_mm;
}bam_vector_pack;

#define PLUS_OR 1
#define MINUS_OR 0


typedef struct barcoded_read_pair{
	unsigned int left : 31;
	unsigned int l_or : 1;
	unsigned int right :31;
	unsigned int r_or :1;
	unsigned long barcode;
	int l_chr :16;
	int r_chr :16;
}barcoded_read_pair;

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
        short **depths;
	double *depth_mean;
	double *depth_std;
	int read_count;
	bit_set_t *chro_bs;
} bam_info;

bam_info *get_bam_info(sonic *);
bam_stats *calculate_bam_statistics(bam_info *, char *bam_path, int number_of_reads);
bam_vector_pack *make_bam_vector_pack();
int read_pcs_bam( bam_info* in_bam, char *bam_path, int pool_no, bam_vector_pack *pack);


bam_vector_pack *read_10X_chr_intra( bam_info* in_bam, char *bam_path, sonic *snc,int chr, bam_stats *);
bam_vector_pack *read_10X_chr_inter( bam_info* in_bam, char *bam_path, sonic *snc,int chr, bam_stats *);
bam_vector_pack *read_10X_chr( bam_info* in_bam, char *bam_path, sonic *snc,int chr, bam_stats *);
bam_vector_pack **read_10X_bam( bam_info* in_bam, char *bam_path, sonic *snc);
void destroy_bams( bam_vector_pack* reads);
#endif

