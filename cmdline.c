#include "cmdline.h"
#include "common.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vector.h"
#include "valorconfig.h"
#include "set.h"

double SV_OVERLAP_RATIO;
int MAX_FRAG_SIZE;
int VALOR_FILTER_GAP;
int VALOR_FILTER_SAT;
int BARCODE_LEN;
double CLONE_MEAN;
double CLONE_STD_DEV;
int MOLECULE_EXT;
long CLONE_MAX_DIST;
int CLONE_MIN_DIST;
int MOLECULE_BIN_SIZE;
int INV_MIN_SIZE;
int INV_MAX_SIZE;
int INVERSION_MIN_REQUIRED_SUPPORT;
int INVERSION_MIN_CLUSTER_SIZE;
int DUP_MIN_SIZE;
int DUP_MAX_SIZE;
int DUP_MAX_DIST;
int DUP_MIN_DIST;
int TRA_MIN_SIZE;
int TRA_MAX_SIZE;
char* VALOR_MOBILE_ELEMENTS;
int DUPLICATION_MIN_CLUSTER_SIZE;
int TRANSLOCATION_MIN_CLUSTER_SIZE;
int DUPLICATION_MIN_REQUIRED_SUPPORT;
int TANDEM_DUPLICATION_MIN_CLUSTER_SIZE;
int TANDEM_DUPLICATION_MIN_SUPPORT;
int DELETION_MIN_REQUIRED_SUPPORT;
int DELETION_MIN_CLUSTER_SIZE;
int MIN_INTER_CLUSTER_SIZE;
int TRA_MIN_INTRA_SPLIT;
double QCLIQUE_LAMBDA;
double QCLIQUE_GAMMA;

int MIN_COVERAGE;
int MAX_COVERAGE;
int EXTENSION;
int MIN_REQUIRED_READS_IN_MOLECULE;
int MIN_QUAL;
int MAX_ALTERNATIVE_CHECK_QUAL;
int READ_SAMPLE_SIZE;
int FILTER1XK;
int MAX_SUPPORT;
int ALTERNATIVE_MAPPING_BIT;
char* ALTERNATIVE_MAPPING_FLAG;
int CHECK_ALTERNATIVE_MAPPINGS;

parameters *get_params(void){
	static parameters *params = NULL;
	if(params == NULL){ params = malloc(sizeof(parameters));}
	return params;
}
void free_params(void /*parameters*/ *vp){
	parameters *p = vp;
	free(p->sonic_file);
	free(p->logfile);
	free(p->outprefix);
	free(p->bam_file);
	free(p);
}
int vstrncmp(const void *v1, const void *v2, size_t len){
    return strncmp((const char*) v1, (const char*) v2, len);
}

sv_type parse_svs( const char * optt){

	vector_t *sv_strings = dang_string_tokenize(optt, ",");
	sv_type sv_to_find = 0;

	int i;
	for(i=0;i<sv_strings->size;i++){
		sv_to_find|=atosv(vector_get(sv_strings,i));
	}
	vector_free(sv_strings);
	return sv_to_find;
}

// #### TODO move these to params ####
void set_valor_option(parameters *param, const char *oc, const char *optarg){
    if(!strcmp(oc,"input")){
        param->bam_file = malloc(strlen(optarg)+1);
        strcpy(param->bam_file, optarg);
    }
    else if(!strcmp(oc,"sonic")){
        param->sonic_file = malloc(strlen(optarg)+1);
        strcpy(param->sonic_file, optarg);
    }
    else if(!strcmp(oc,"out")){
        param->outprefix = malloc(strlen(optarg)+1);
        strcpy(param->outprefix, optarg);
    }
    else if(!strcmp(oc,"minimum-split-distance")){
        CLONE_MIN_DIST = atoi(optarg);
    }
    else if(!strcmp(oc,"maximum-split-distance")){
        CLONE_MAX_DIST = atol(optarg);
    }
    else if(!strcmp(oc,"molecule-bin-size")){
        MOLECULE_BIN_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-min-size")){
        DUP_MIN_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-max-size")){
        DUP_MAX_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-min-dist")){
        DUP_MIN_DIST = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-max-dist")){
        DUP_MAX_DIST = atoi(optarg);
    }
    else if(!strcmp(oc,"deletion-min-support")){
        DELETION_MIN_REQUIRED_SUPPORT = atoi(optarg);
    }
    else if(!strcmp(oc,"deletion-min-cluster-size")){
        DELETION_MIN_CLUSTER_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-min-support")){
        DUPLICATION_MIN_REQUIRED_SUPPORT = atoi(optarg);
    }
    else if(!strcmp(oc,"duplication-min-cluster-size")){
        DUPLICATION_MIN_CLUSTER_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"inversion-min-size")){
        INV_MIN_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"inversion-max-size")){
        INV_MAX_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"inversion-min-support")){
        INVERSION_MIN_REQUIRED_SUPPORT = atoi(optarg);
    }
    else if(!strcmp(oc,"translocation-min-cluster-size")){
        TRANSLOCATION_MIN_CLUSTER_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"inversion-min-cluster-size")){
        INVERSION_MIN_CLUSTER_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"translocation-min-size")){
        TRA_MIN_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"translocation-max-size")){
        TRA_MAX_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc,"sv-overlap-ratio")){
        SV_OVERLAP_RATIO = atof (optarg);
    }
    else if(!strcmp(oc,"mobile-elements")){
        VALOR_MOBILE_ELEMENTS = malloc(strlen(optarg) + 1);
        strcpy(VALOR_MOBILE_ELEMENTS, optarg);
    }
    else if(!strcmp(oc, "no-gap-filter")){
        VALOR_FILTER_GAP = 0;
    }
    else if(!strcmp(oc, "no-satellite-filter")){
        VALOR_FILTER_SAT = 0;
    }
    else if(!strcmp(oc, "max-frag-size")){
        MAX_FRAG_SIZE = atoi(optarg);
    }
    else if(!strcmp(oc, "barcode-len")){
        BARCODE_LEN = atoi(optarg);
    }
    else if(!strcmp(oc, "minimum-molecule-coverage")){
        MIN_COVERAGE = atoi(optarg);
    }
    else if(!strcmp(oc, "maximum-molecule-coverage")){
        MAX_COVERAGE = atoi(optarg);
    }
    else if(!strcmp(oc, "quasi-clique-lambda")){
        QCLIQUE_LAMBDA = atof(optarg);
    }
    else if(!strcmp(oc, "quasi-clique-gamma")){
        QCLIQUE_GAMMA = atof(optarg);
    }
    else if(!strcmp(oc, "molecule-extension-distance")){
        EXTENSION = atoi(optarg);
    }
    else if(!strcmp(oc, "min-required-reads-in-molecule")){
        MIN_REQUIRED_READS_IN_MOLECULE = atoi(optarg);
    }
    else if(!strcmp(oc, "min-alignment-quality")){
        MIN_QUAL = atoi(optarg);
    }
    else if(!strcmp(oc, "single-copy-chr")){
        // Handle outside
    }
    else if(!strcmp(oc, "svs-to-find")){
        param->svs_to_find = parse_svs(optarg);
    }
    else if(!strcmp(oc, "ploidy")){
        param->ploidy = atoi(optarg);
    }
    else if(!strcmp(oc, "log-file")){
        param->logfile = malloc(strlen(optarg)+1);
        strcpy(param->logfile, optarg);
    }
    else if(!strcmp(oc, "contig-count")){
        param->chromosome_count = atoi(optarg); 
    }
    else if(!strcmp(oc, "threads")){
        param->threads = atoi(optarg);
    }
    else if(!strcmp(oc, "preset")){
        if(strcmp(optarg, "None")){
            fprintf(stderr, "Presets are not yet implemented!\n");
            exit(-1);
        }
    }
    else if(!strcmp(oc, "read-coverage")){
        param->coverage = atoi(optarg);
    }
    else if(!strcmp(oc,"molecule-discovery-window" )){
        MOLECULE_EXT = atoi(optarg);
    }
    else if(!strcmp(oc,"stats-sample-size")){
       READ_SAMPLE_SIZE = atoi(optarg); 
    }
    else if(!strcmp(oc,"max-discordant-support")){
       MAX_SUPPORT = atoi(optarg); 
    }
    else{ //This should only run when setting the defaults, if the option is not handled. getopt will not let invalid options to be evaluated.
        fprintf(stderr, "Unknown argument %s with value %s. Please revise the running command!\n", oc, optarg);
        exit(-1);
    }
}





typedef struct {
    char *name;
    int         has_arg;
    int        *flag;
    int         val;
} optclone;

typedef struct opt_and_help{
    optclone o;
    char *default_value;
    char *help;
} opt_and_help;


void free_opt_and_help( void *voah){
    opt_and_help *oah = voah;
    free(oah->o.name);
    free(oah->default_value);
    free(oah->help);
    free(oah);
}

typedef struct arg_manager{
    vector_t *mandatory;
    vector_t *optional;
    vector_t *opt_str_builder;
    char *header;
    char *footer;
} arg_manager;

arg_manager *make_argument_manager(char *header, char *footer){
    arg_manager *argm = malloc(sizeof(arg_manager));
    argm->mandatory = vector_init(sizeof(opt_and_help), 16);
    argm->optional = vector_init(sizeof(opt_and_help), 16);
    

    vector_set_remove_function(argm->mandatory, free_opt_and_help);
    vector_set_remove_function(argm->optional, free_opt_and_help);

    argm->opt_str_builder = vector_init(sizeof(char), 16);

    argm->header = malloc(strlen(header)+1);
    argm->footer = malloc(strlen(footer)+1);
    strcpy(argm->header, header);
    strcpy(argm->footer, footer);
    return argm;
}

void arg_manager_set_defeaults(arg_manager *argm, parameters * param){
    int i;
    for(i = 0; i < argm->optional->size; ++i){
        opt_and_help *oah = vector_get(argm->optional, i);
        if( oah->default_value != 0 && strcmp(oah->o.name,"") != 0 && strcmp(oah->default_value, "") != 0){
            set_valor_option(param, oah->o.name, oah->default_value);
        }
    }
}

void free_argument_manager(arg_manager *argm){
    vector_free(argm->mandatory);
    vector_free(argm->optional);
    vector_free(argm->opt_str_builder);
    free(argm->header);
    free(argm->footer);
    free(argm);
}

#define MANDATORY_ARG 0

void add_short_argument( arg_manager *argm, int has_arg, int sarg, char *help, char* default_value){

    vector_t *avec;
    opt_and_help *optt = malloc(sizeof(opt_and_help));

    optt->o.name    = malloc(1);
    optt->o.name[0] = 0;
    optt->o.has_arg = has_arg;
    optt->o.flag    = 0;
    optt->o.val     = sarg;
    optt->default_value = 0;
  

    char colon = ':';
    vector_put(argm->opt_str_builder, &sarg);
    if(has_arg == required_argument){
        vector_put(argm->opt_str_builder, &colon);
    }
    else if(has_arg == optional_argument){
        vector_put(argm->opt_str_builder, &colon);
        vector_put(argm->opt_str_builder, &colon);
    }
    if(default_value == 0){
        avec = argm->mandatory;
    }else{
        avec = argm->optional;
        optt->default_value = malloc(strlen(default_value)+1);
        strcpy(optt->default_value, default_value);
    }
    
    optt->help       = malloc(strlen(help)+1);
    strcpy(optt->help, help);

    vector_soft_put(avec, optt);

}

void add_long_argument( arg_manager *argm, const char *larg, int has_arg, int *flag, int sarg, char *help, char* default_value){
    vector_t *avec;
    opt_and_help *optt = malloc(sizeof(opt_and_help));

    int larg_size = strlen(larg);
    optt->o.name    = malloc(larg_size + 1);
    strncpy(optt->o.name, larg, larg_size);
    optt->o.name[larg_size]=0; 
    optt->o.has_arg = has_arg;
    optt->o.flag    = flag;
    optt->o.val     = sarg;
    optt->default_value = 0;

    if( sarg <= 'z'){
        char colon = ':';
        vector_put(argm->opt_str_builder, &sarg);
        if(has_arg == required_argument){
            vector_put(argm->opt_str_builder, &colon);
        }
        else if(has_arg == optional_argument){
            vector_put(argm->opt_str_builder, &colon);
            vector_put(argm->opt_str_builder, &colon);
        }
    }
    if(default_value == 0){
        avec = argm->mandatory;
    }else{
        avec = argm->optional;
        optt->default_value = malloc(strlen(default_value)+1);
        strcpy(optt->default_value, default_value);
    }
    optt->help       = malloc(strlen(help)+1);
    strcpy(optt->help, help);

    vector_soft_put(avec, optt);
}

void print_help(FILE *stream, arg_manager *argm){
    fprintf(stream,  "%s\n",argm->header);
    fprintf(stream,  "====== Required Arguments ======\n\n");
    int i;
    for(i = 0; i < argm->mandatory->size; ++i){
        opt_and_help *oah = vector_get(argm->mandatory, i);
        if( oah->o.val <= 122 && oah->o.val >=64){
            fprintf(stream, "  -%c,\t",oah->o.val);
        }
        else{
            fprintf(stream, "\t ");
        }
        if( strncmp(oah->o.name, "",256)){
            fprintf(stream, "--%s,",oah->o.name);
        }else{
            fprintf(stream, "\t\t");
        }
        fprintf(stream, "\t%s\n", oah->help);
    }

    fprintf(stream,  "\n=========== Options ============\n\n");
    for(i = 0; i < argm->optional->size; ++i){
        opt_and_help *oah = vector_get(argm->optional, i);
        if( oah->o.val <= 122 && oah->o.val >=64){
            fprintf(stream, "  -%c,\t",oah->o.val);
        }
        else{
            fprintf(stream, "\t ");
        }
        if( strncmp(oah->o.name, "",256)){
            fprintf(stream, "--%s,",oah->o.name);
        }else{
            fprintf(stream, "\t\t");
        }
        fprintf(stream, "\t%s", oah->help);
        if(oah->default_value[0] != 0){
            fprintf(stream, " [%s]\n", oah->default_value);
        }else{
            fprintf(stream, "\n");
        }
    }
    fprintf(stream, "\n================================\n");
    fprintf(stream,  "%s\n",argm->footer);
}


parameters *parse_args(int argc, char **argv){
    parameters *param = get_params();
    memset(param,0, sizeof(parameters));
    char *valor_art =   "\t\t┌────────────────────────┐ \n"
                        "\t\t│         VALOR2         │▒\n"
                        "\t\t╞════════════════════════╡▒\n"
                        "\t\t│    VAriation with      │▒\n"
                        "\t\t│ LOng Range information │▒\n"
                        "\t\t└────────────────────────┘▒\n"
                        "\t\t ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒\n";
    
    arg_manager * argm = make_argument_manager(valor_art, "Send your issues at https://github.com/BilkentCompGen/valor/issues");
    add_long_argument( argm, "input", required_argument, 0, 'i', "Input sorted linked-read bam file", MANDATORY_ARG);
    add_long_argument( argm, "sonic", required_argument, 0, 's', "Sonic file build from the genome reference", MANDATORY_ARG);
    add_long_argument( argm, "out", required_argument, 0, 'o', "Output prefix", MANDATORY_ARG);
        
    add_long_argument( argm, "svs-to-find", required_argument, 0, 'f', "Comma separated list of SVs to find (or ALL). e.g. DUP,INV,DEL", "ALL");
    add_long_argument( argm, "ploidy", required_argument, 0, 'p', "Number of copies per chromosome", "2");
    add_long_argument( argm, "single-copy-chr", required_argument, 0, 'y', "Comma separated list of single copy chromosomes", "");
    add_long_argument( argm, "log-file", required_argument, 0, 'l', "Log file path", "valor.log");
    add_long_argument( argm, "contig-count", required_argument, 0, 'c', "First N contigs to run VALOR2", "24");

    add_long_argument( argm, "threads", required_argument, 0, 't', "Number of threads", "8");
    add_long_argument( argm, "preset" , required_argument, 0, 'x', "Presets" , "None");



    add_long_argument( argm, "barcode-len" , required_argument, 0, 'b', "barcode length" , "16");
    add_long_argument( argm, "read-coverage", required_argument, 0, 'z', "Expected read coverage of the input bam file", "40");

    int ai = 123;

    add_long_argument( argm, "no-gap-filter", no_argument, 0, ai++, "Dont Filter Gaps (This will slow down VALOR2 significantly, don't use if you don't know what you are doing.)", "False");
    add_long_argument( argm, "no-satellite-filter", no_argument, 0, ai++, "Dont Filter Satellites (This will slow down VALOR2 significantly, don't use if you don't know what you are doing.)", "False");
    add_long_argument( argm, "molecule-discovery-window" , required_argument, 0, ai++, "Initial molecule discovery window" , "80000");
    add_long_argument( argm, "molecule-extension-distance" , required_argument, 0, ai++, "Extension distance after initial discovery window" , "10000");
    add_long_argument( argm, "minimum-split-distance" , required_argument, 0, ai++, "Minimum split distance (2 * Expected Molecule Size)" , "80000");
    add_long_argument( argm, "maximum-split-distance" , required_argument, 0, ai++, "Maximum split distance" , "3000000000");
    add_long_argument( argm, "molecule-bin-size" , required_argument, 0, ai++, "Molecule bin size (For molecule coverage calculation)" , "8000");

    add_long_argument( argm, "duplication-min-size" , required_argument, 0, ai++, "Minimum duplication size to discover " , "80000");
    add_long_argument( argm, "duplication-max-size" , required_argument, 0, ai++, "Maximum duplication size to discover " , "15000000");
    add_long_argument( argm, "duplication-min-dist" , required_argument, 0, ai++, "Minimum duplication distance to discover " , "80000");
    add_long_argument( argm, "duplication-max-dist" , required_argument, 0, ai++, "Maximum duplication distance to discover (Smaller == Faster) " , "15000000");
    add_long_argument( argm, "duplication-min-support" , required_argument, 0, ai++, "Minimum short read support for an duplication event" , "24");
    add_long_argument( argm, "duplication-min-cluster-size" , required_argument, 0, ai++, "Minimum number of unique duplication split molecule pairs" , "24");

    add_long_argument( argm, "inversion-min-size" , required_argument, 0, ai++, "Minimum inversion size to discover " , "80000");
    add_long_argument( argm, "inversion-max-size" , required_argument, 0, ai++, "Maximum inversion size to discover (Smaller == Faster) " , "15000000");

    add_long_argument( argm, "inversion-min-support" , required_argument, 0, ai++, "Minimum short read support for an inversion event" , "12");
    add_long_argument( argm, "inversion-min-cluster-size" , required_argument, 0, ai++, "Minimum number of unique inversion split molecule pairs" , "24");

    add_long_argument( argm, "deletion-min-support" , required_argument, 0, ai++, "Minimum short read support for an deletion event" , "8");
    add_long_argument( argm, "deletion-min-cluster-size" , required_argument, 0, ai++, "Minimum number of unique deletion split molecules" , "8");
   
    add_long_argument( argm, "translocation-min-size" , required_argument, 0, ai++, "Minimum translocation size to discover " , "80000");
    add_long_argument( argm, "translocation-max-size" , required_argument, 0, ai++, "Maximum translocation size to discover (Smaller == Faster) " , "15000000");
    add_long_argument( argm, "translocation-min-cluster-size" , required_argument, 0, ai++, "Minimum number of unique translocation split molecules" , "16");
    
    add_long_argument( argm, "sv-overlap-ratio" , required_argument, 0, ai++, "Minimum molecule overlap ratio to consider 2 SVs overlapping" , "0.25");

    add_long_argument( argm, "mobile-elements" , required_argument, 0, ai++, "colon seperated list of mobile elements to consider" , "Alu:L1:SVA:HERV");
    
    add_long_argument( argm, "quasi-clique-lambda" , required_argument, 0, ai++, "Lambda value of quasi clique approximation between 0 and 1" , "0.5");
    add_long_argument( argm, "quasi-clique-gamma"  , required_argument, 0, ai++, "Gamma  value of quasi clique approximation between 0 and 1" , "0.6");
    
    add_long_argument( argm, "minimum-molecule-coverage"  , required_argument, 0, ai++, "Minimum coverage for a valid molecule", "0");

    add_long_argument( argm, "maximum-molecule-coverage"  , required_argument, 0, ai++, "Maximum coverage for a valid molecule", "50");
    add_long_argument( argm, "min-required-reads-in-molecule" , required_argument, 0, ai++, "Minimum number of reads per molecule", "7");

    add_long_argument( argm, "max-frag-size"  , required_argument, 0, ai++, "Maximum read fragment size", "1000");
    add_long_argument( argm, "max-discordant-support"  , required_argument, 0, ai++, "Maximum read fragment size", "100");


    add_long_argument( argm, "min-alignment-quality", required_argument, 0, ai++, "Minimum alignment quality to a read to be used", "1");
    add_long_argument( argm, "stats-sample-size", required_argument, 0, ai++, "Number of read pairs to calculate read statistics", "1000000");

    add_short_argument( argm, no_argument, 'v', "Verbose", "");
    add_long_argument( argm, "help", no_argument, 0, 'h', "Shows this help screen", "");
    add_long_argument( argm, "version", no_argument, 0, 'V', "Print VALOR2 version", "");

    arg_manager_set_defeaults(argm, param);

    size_t num_options = argm->mandatory->size + argm->optional->size;
    struct option *long_options = malloc(sizeof(struct option) * (num_options + 1)); //+1 is required by null option
    int i;
    int j;
    
    if( argc == 1){
        print_help(stderr,argm);
        free(long_options);
        free_argument_manager(argm);

        exit(-1);
    }
    char *opt_string = malloc(argm->opt_str_builder->size + 1);
    for(i = 0; i< argm->opt_str_builder->size; ++i){
        opt_string[i] = *(char *) vector_get(argm->opt_str_builder,i);
    }

    opt_string[i] = 0;

    for(i = 0, j = 0; i < argm->mandatory->size; i+=1){
        struct option *ooo =vector_get(argm->mandatory,i);
        if(ooo->name == 0){ continue;};
        long_options[j] = *ooo;
        ++j;
    }
    for(i = 0; i < argm->optional->size; i+=1){
        struct option *ooo =vector_get(argm->optional,i);
        if(ooo->name == 0){ continue;};
        long_options[j] = *ooo;
        ++j;
    }

    long_options[j] = (struct option){0,0,0,0};

    int verbosity = 0;
    int o;
    int index = 0;
    param->haplotype_chrs = vector_init(sizeof(char)*128,128);
    char *buffer;
    const char *oc = 0;

    set_t *used_args = set_init(128,0);
    used_args->key_cmp = vstrncmp;
    used_args->hf = SuperFastStringHash;
    //printf("%s\n",opt_string);
//// i:s:o:f:p:y:l:c:t:x:gvb:z:{:|:}:~::▒:▒:▒:▒:▒:▒:▒:▒:▒:▒:▒:▒:▒:hV
    while( ( o = getopt_long( argc, argv, opt_string, long_options, &index)) != -1){
        switch( o){
            case 'y':
                buffer = malloc(128);
				strncpy(buffer,optarg,127);
				vector_soft_put(param->haplotype_chrs,buffer);
            break;
            case 'h':
                print_help(stdout,argm);
                free(long_options);
                free_argument_manager(argm);
                free(opt_string);
                exit(1);
            case 'v':
                ++verbosity;
            break;
            case '?':
                exit(1);
            break;
            default:
                if( index >= num_options || index < 0){
                    fprintf(stderr, "Wrong index : %d, with value %s, with o: %d\n",index,optarg,o);
                }
                for(i=0;i<num_options;++i){
                    if(long_options[i].val == o){
                        index = i;
                        break;
                    }

                }
   
                oc = long_options[index].name;
                char *arg_name = malloc(strlen(oc) + 1);
				strncpy(arg_name,oc,strlen(oc));
                arg_name[strlen(oc)]=0;
                printf("%d, %s: %s\n",o,arg_name,optarg);
                set_soft_put(used_args,arg_name);
                //printf("%lu\n",used_args->number_of_items);
                set_valor_option(param, oc, optarg);
            break;
        }

    }
    free(opt_string);
    int failed = 0;
    for(i = 0; i < argm->mandatory->size; i+=1){
        struct option *ooo =vector_get(argm->mandatory,i);
        if(ooo->name == 0){ continue;};
        if(!set_has(used_args,ooo->name)){
            fprintf(stderr,"\n%s is a mandatory argument!!!\n",ooo->name);
            failed = 1;
        }

    }
    //FREE
    set_free(used_args);
    if( failed){
        print_help(stderr,argm);

    }    
    free(long_options);
    free_argument_manager(argm);
    if( failed){
        exit(-1);
    }    
    return param;
}

