#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#include "common.h"
#include "cmdline.h"
#include "config.h"
#include "valor.h"
#include "readbam.h"
#include "readbed.h"
#include "statistics.h"
#include "recovermolecules.h"
#include "interval10X.h"
#include "structural_variation.h"
#include "interc_sv.h"
#include "vector.h"
#include "bitset.h"
#include "clique.h"
#include "cluster.h"
#include "cnv.h"
#include "sonic/sonic.h"
#include "graph.h"
#include "progress.h"
#include <omp.h>

FILE *logFile = NULL;
double CLONE_MEAN;
double CLONE_STD_DEV;

int main( int argc, char **argv){

    parameters * params = init_params();
    if(parse_command_line(argc,argv,params)){
        return 0;
    }
#ifdef _OPENMP
    omp_set_num_threads(params->threads);
#endif
    sv_type svs_to_find =  params->svs_to_find;

    char *bamname = params->bam_file;;

    int i,j,k;
    time_t rawtime;
    struct tm *timeinfo;


    char *out_file_path = malloc((strlen("-predicted_svs.bedpe")+strlen(params->outprefix)+1)*sizeof(char));
    sprintf(out_file_path,"%s-predicted_svs.bedpe",params->outprefix);
    FILE *outbedfile = fopen(out_file_path,"w+");
    free(out_file_path);
    time( &rawtime);

    timeinfo = localtime( &rawtime);


    printf("\n\nVALOR: Variation with LOng Range\n");
    printf("Version: %s\n", VALOR_VERSION);
    printf("Build Date: %s\n",BUILD_DATE);
    printf("Output: %s-predicted_svs.bedpe\n",params->outprefix);
    printf("Logfile: %s\n", params->logfile);


    sonic *snc = sonic_load(params->sonic_file);
    printf("Reading BAM file: %s\n", bamname);

    logFile = safe_fopen(params->logfile,"w+");                                                                    
    printvalorconfig(logFile);
    char *molecule_bed_path = malloc((strlen(params->outprefix) + strlen("-molecules.bed") + 1) * sizeof(char));
    sprintf(molecule_bed_path,"%s-molecules.bed",params->outprefix);                                            
    FILE *reset_molecule_bed = fopen(molecule_bed_path,"w+");
    fclose(reset_molecule_bed);
    //////                                                                                                        
    //
    bam_info *in_bams = get_bam_info(snc);

    fprintf( logFile, "#CreationDate=%d.%d.%d\n\n",
            timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday); 
    in_bams->sample_name = NULL;


    bam_stats *stats = calculate_bam_statistics(in_bams, bamname, READ_SAMPLE_SIZE);

    //	bam_vector_pack **reads = read_10X_bam(in_bams,bamname,snc);
    bam_vector_pack **reads = malloc(sizeof(bam_vector_pack) * snc->number_of_chromosomes);//read_10X_bam(in_bams,bamname,snc);


    int first_skipped = 0;
    for( i = 0; i < params->chromosome_count ;i++){


        reads[i] = read_10X_chr_intra(in_bams,bamname,snc,i,stats);
        if(reads[i]->concordants->size == 0){
            destroy_intra_bams(reads[i]);
            printf("No reads for chromosome %s %s %s.\r",
                    snc->chromosome_names[first_skipped],
                    (i-first_skipped==1?"and":"to"),
                    snc->chromosome_names[i]);
            continue;
        }
    



        bit_set_set_bit(get_bam_info(NULL)->chro_bs,i,1);
        printf("\nFinding structural variants in chromosome %s\n",snc->chromosome_names[i]);
        first_skipped = i+1;

        printf("Recovering split molecules..\n");
        vector_t *regions = recover_molecules(reads[i]->concordants);

        if(params->svs_to_find & SV_TRANSLOCATION){
            append_molecules_to_bed(regions,molecule_bed_path,i);
        }
        in_bams->depths[i] = make_molecule_depth_array(regions,snc,i);

        in_bams->depth_mean[i] = make_global_molecule_mean(in_bams->depths[i],snc,i);
        in_bams->depth_std[i] = make_global_molecule_std_dev(in_bams->depths[i],snc,i,in_bams->depth_mean[i]);
        in_bams->depth_std[i] = MIN(in_bams->depth_std[i],in_bams->depth_mean[i]/2);
        printf("Global molecule depth mean: %lf\nGlobal molecule depth standard deviation: %lf\n",
                in_bams->depth_mean[i],in_bams->depth_std[i]);

        VALOR_LOG("chr: %s\nmolecule mean depth: %lf\nmolecule std-dev depth: %lf\n", snc->chromosome_names[i], in_bams->depth_mean[i], in_bams->depth_std[i]);
        VALOR_LOG("Initial molecule count: %zu\n",regions->size);
        filter_molecules(regions,snc,i);

        VALOR_LOG("Filtered molecule count: %zu\n",regions->size);
        qsort(regions->items,regions->size,sizeof(void*),interval_start_comp);  
        vector_t *groups = group_overlapping_molecules(regions);
        groups->rmv = &vector_free;

        CLONE_MEAN = molecule_mean(regions);
        CLONE_STD_DEV = molecule_group_std(groups,CLONE_MEAN);
        VALOR_LOG("Molecule size mean: %lf\nMolecule size std-dev: %lf\n", CLONE_MEAN, CLONE_STD_DEV);
        vector_free(groups);

        qsort(regions->items,regions->size,sizeof(void*),barcode_comp);  

        printf("Molecule Count: %zu\tMolecule mean: %lf\tMolecule std-dev: %lf\n",regions->size,CLONE_MEAN,CLONE_STD_DEV);

        vector_t *split_molecules = discover_split_molecules(regions);
        vector_free(regions);

        VALOR_LOG("Split molecule candidate count: %zu\n",split_molecules->size);

        printf("Matching split molecules\n");
        vector_t *variations = find_svs(split_molecules,svs_to_find, i);

        VALOR_LOG("Matched split molecule pair count: %zu\n", variations->size);

        qsort(split_molecules->items,split_molecules->size,sizeof(void*),interval_pair_comp);  

        printf("%zu candidate variations are made\n",variations->size);
        update_sv_supports_b(variations,
                reads[i]);
        variations->REMOVE_POLICY = REMP_LAZY;

        filter_unsupported_pm_splits(split_molecules, reads[i]->pm_discordants);
        vector_filter(variations,sv_is_proper);

        VALOR_LOG("Structural variation candidate count: %zu\n",variations->size);

        printf("%zu candidate variations are left after filtering\n",variations->size);
#if FILTER1XK //If number of variations are more than MAX_INVERSIONS_IN_GRAPH, do random selection
        if(variations->size > MAX_INVERSIONS_IN_GRAPH){
            srand(time(NULL));
            variations->REMOVE_POLICY=REMP_LAZY;
            for(k=0;k<variations->size;k++){
                if(rand()%variations->size>MAX_INVERSIONS_IN_GRAPH){
                    vector_remove(variations,k);
                }
            }
            vector_defragment(variations);

            printf("%zu candidate variations are left after random selection to decrease size\n",variations->size);
        }
#endif

        variations->REMOVE_POLICY=REMP_SORTED;

        graph_t *sv_graph = make_sv_graph(variations);
        /*
           fff = fopen("sv_graph_out.out","w+");
           graph_print(sv_graph,fff);
           fclose(fff);

*/
        printf("Finding SV Clusters\n\n");
        vector_t *clusters = vector_init(sizeof(sv_cluster),50);
        vector_set_remove_function(clusters,&sv_cluster_destroy);

        graph_trim(sv_graph);

        vector_t *comps = sv_g_dfs_components(sv_graph);

        if(comps->size == 0){continue;}

        for(k=0;k<comps->size;k++){
            vector_t *garbage = vector_get(comps,k);
            size_t initial_size = garbage->size;
            qsort(garbage->items, garbage->size, sizeof(sv_t *),&sv_comp);
            sv_type _type = ((sv_t *)vector_get(garbage,0))->type;
            if(garbage->size < what_is_min_cluster_size(_type,params->chr_copy_count[i])){continue;}
            graph_t *garbage_graph = sv_graph;//= make_sv_graph(garbage);
            int iteration_no = 0;
            while(garbage_graph->number_of_items > 2){
                if(garbage->size < initial_size /2 ){break;}

                clique_t *c = clique_find_clique(garbage_graph,garbage,0,params->quasi_clique_lambda,params->quasi_clique_gamma);

                if(c==NULL||c->v_prime<=0){clique_free(c);break;}
                sv_cluster *svc_garbage = sv_cluster_make(c);
                clique_free(c);

                if(svc_garbage==NULL){break;}
                if(svc_garbage->items->size < what_is_min_cluster_size(_type,params->chr_copy_count[i])){
                    sv_graph_reset(garbage_graph);
                    sv_cluster_graph_fix(svc_garbage,garbage,sv_graph);
                    sv_cluster_destroy(svc_garbage);
                    continue;
                }


                sv_graph_reset(garbage_graph);

                sv_cluster_graph_fix(svc_garbage,garbage,garbage_graph);



                vector_soft_put(clusters,svc_garbage);
                iteration_no++;
            }
            //			vector_free(garbage);
            //	graph_free(garbage_graph);
        }

        printf("Clustering is finished, found %zu variant clusters\n",clusters->size);

        vector_free(comps);

        qsort(clusters->items, clusters->size, sizeof( void*), cluster_comp);
        printf("Printing variant calls\n");
        for(j=0;j<clusters->size;j++){
            sv_cluster *svc = vector_get(clusters,j);


            sv_t *first = vector_get(svc->items,0);
            if(svc->items->size < what_is_min_cluster_size(first->type,params->chr_copy_count[i])){
                fprintf(logFile,"%s\t%d\t%d\t%s\t%d\t%d\t%s\t%zu\t%d\tSmall CLSTR\n",
                    snc->chromosome_names[i],
                    svc->break_points->start1,
                    svc->break_points->end1,
                    snc->chromosome_names[i],
                    svc->break_points->start2,
                    svc->break_points->end2,
                    sv_type_name(first->type),
                    svc->items->size,
                    svc->supports[0]+svc->supports[1]
                   );

                
                continue;
            }


            double mean_depth = 0;
            if(first->type == SV_DIRECT_DUPLICATION || first->type == SV_INVERTED_DUPLICATION){
                if(first->orientation == DUP_FORW_COPY){
                    mean_depth = get_depth_region(in_bams->depths[i],svc->break_points->start1,svc->break_points->end1);
                }else if(first->orientation == DUP_BACK_COPY){
                    mean_depth = get_depth_region(in_bams->depths[i],svc->break_points->start1,svc->break_points->end1);
                }

//                if(mean_depth < 1.5 * in_bams->depth_mean[i]){ continue;}
                //if(mean_depth < in_bams->depth_mean[i] + 1.25 * in_bams->depth_std[i]){ continue;}
            }
            else if (first->type == SV_TRANSLOCATION || first->type == SV_INVERTED_TRANSLOCATION){

                interval_pair deletion_interval;


                deletion_interval = (interval_pair){.start1=svc->break_points->start1-CLONE_MEAN/2,.end1=svc->break_points->start1,.start2=svc->break_points->end1,.end2=svc->break_points->end1+CLONE_MEAN/2,.barcode=0};


                size_t pos = 0;//split_molecule_binary_search(split_molecules,deletion_interval);
                if( pos == -1){ continue;}
                splitmolecule_t *cand = vector_get(split_molecules,pos);

                  vector_t *found_splits = vector_init(sizeof(interval_pair),10);
                while( pos < split_molecules->size){ // && cand->start1 < deletion_interval.end + 50000){

                    if(interval_pair_overlaps(&deletion_interval,cand,CLONE_MEAN/2)){
                         vector_put(found_splits,cand);
                    }

                    /*if( cand->end1 < deletion_interval.start + CLONE_MEAN /2  && 
                            cand->start1 + cand->end1 > 2 * deletion_interval.start - 3* CLONE_MEAN  && 
                            cand->start2 + cand->end2 < 2 * deletion_interval.end +3* CLONE_MEAN &&
                            cand->start2 >deletion_interval.end - CLONE_MEAN /2){

                        flag = 1;
                        break;
                    }*/
                    pos++;
                    cand = vector_get(split_molecules,pos);
                }
                size_t tra_pm_cnt = found_splits->size;
                vector_free(found_splits);
                if(tra_pm_cnt < TRA_MIN_INTRA_SPLIT){
                    vector_free(found_splits);
                    continue;
                }

                mean_depth = get_depth_region(in_bams->depths[i],first->AB.end1, first->AB.start2);
            }else if (first->type == SV_TANDEM_DUPLICATION){
                mean_depth = get_depth_region(in_bams->depths[i],first->AB.start1, first->AB.end2);
            }
            else if (first->type == SV_DELETION){
                mean_depth = get_depth_region(in_bams->depths[i],first->AB.end1, first->AB.start2);
            }else {
                mean_depth = get_depth_region(in_bams->depths[i],first->AB.end1,first->CD.start1)/2 + get_depth_region(in_bams->depths[i],first->AB.end2,first->CD.start2)/2;
            }

            fprintf(outbedfile,"%s\t%d\t%d\t%s\t%d\t%d\t%s\t%zu\t%d\t%lf\n",
                    snc->chromosome_names[i],
                    svc->break_points->start1,
                    svc->break_points->end1,
                    snc->chromosome_names[i],
                    svc->break_points->start2,
                    svc->break_points->end2,
                    sv_type_name(first->type),
                    svc->items->size,
                    svc->supports[0]+svc->supports[1],
                    mean_depth       
                   );
        }

        fflush(outbedfile);

        if((svs_to_find & SV_TRANSLOCATION) == 0){
            destroy_intra_bams(reads[i]);
            free(in_bams->depths[i]);
        }else{
            vector_free(reads[i]->mp_discordants);
            vector_free(reads[i]->mm_discordants);
            vector_free(reads[i]->pp_discordants);
            vector_free(reads[i]->concordants);
            reads[i]->mp_discordants = NULL;
            reads[i]->pp_discordants = NULL;
            reads[i]->mm_discordants = NULL;
            reads[i]->concordants = NULL;
        }

        vector_free(split_molecules);
        vector_free(variations);
        vector_free(clusters);
        graph_free(sv_graph);

        printf("Reading next chromosome\n");
    }
    printf("\n");

    if(svs_to_find & SV_TRANSLOCATION || svs_to_find & SV_INVERTED_TRANSLOCATION){
        printf("Looking for translocations.\n");
        printf("Reading from temp molecule file.\n");
        vector_t **molecules = read_molecules_from_bed(molecule_bed_path);
        vector_t *variants = find_interchromosomal_events_lowmem(molecules,reads, bamname);
        for(k=0;k<variants->size;k++){
            vector_t *sub_vec = vector_get(variants,k);
            for(i=0;i<sub_vec->size;i++){
                inter_sv_call_bed_print(outbedfile,vector_get(sub_vec,i));
            }
        }
        vector_free(variants);
        for(k=0;k<snc->number_of_chromosomes;k++){
            vector_free(molecules[k]);
        }
        free(molecules);	
 
    }
    
    free(molecule_bed_path);
    free(in_bams->depths);
    free(in_bams->depth_mean);
    free(in_bams->depth_std);
    bit_set_free(in_bams->chro_bs);
    freeMem(in_bams, sizeof(bam_info));
    free(stats);



    free(reads);
    free_sonic(snc);
    fclose(logFile);
    fclose(outbedfile);
    free_params(params);
    return 0;
}
