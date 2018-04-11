#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


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
#include "vector.h"
#include "clique.h"
#include "cluster.h"
#include "cnv.h"
#include "sonic/sonic.h"
#include "graph.h"
#include "progress.h"
int CUR_CHR = -1;
FILE *logFile = NULL;
double CLONE_MEAN;
double CLONE_STD_DEV;



int main( int argc, char **argv){

	sv_type svs_to_find = SV_INVERTED_DUPLICATION;
	bam_info *in_bams;
	char *bamname = argv[1];

	int i,j,k;
	time_t rawtime;
	struct tm *timeinfo;
	vector_t **regions; // Vector of Inverval_10X
	vector_t *groups;   // Vector of Interval_10X
	vector_t **variations; //Vector of sv_t
	vector_t **clusters;   //Vector of cluster_t

	//TODO Take sonic as a parameter

	//	sonic *snc = sonic_load("aux/ucsc_hg19.sonic");
	sonic *snc = sonic_load(argv[2]);
	FILE *outbedfile = fopen(OUT_DIR"/predicted_svs.bed","w+");

	time( &rawtime);

	timeinfo = localtime( &rawtime);

	printf("Output Directory: "OUT_DIR"\n");
	printf("Logfile name: "VALOR_LOG_FILE"\n");
	print_quote();

	logFile = safe_fopen(VALOR_LOG_FILE,"w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n",
			timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday); 

	in_bams = (bam_info *) getMem( sizeof(bam_info));
	in_bams->sample_name = NULL;


	bam_stats *stats = calculate_bam_statistics(in_bams, bamname, READ_SAMPLE_SIZE);

//	bam_vector_pack **reads = read_10X_bam(in_bams,bamname,snc);
	bam_vector_pack **reads = malloc(sizeof(bam_vector_pack) * snc->number_of_chromosomes);//read_10X_bam(in_bams,bamname,snc);

	regions = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	variations = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	clusters = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	int first_skipped = 0;
	for( i = 0; i < 24 /*snc->number_of_chromosomes*/;i++){
#if DEVELOPMENT_
		char fileculname[255];
		sprintf(fileculname,OUT_DIR"clusters%d.out",i+1);
		FILE *filecul = fopen( fileculname,"w+");
		if(filecul==NULL){ fprintf(stderr,"Could not open %s\n",fileculname); return -1;}
#endif
		CUR_CHR = i;
		reads[i] = read_10X_chr(in_bams,bamname,snc,i,stats);
		if(reads[i]->concordants->size == 0){
			printf("No Reads for Chromosome %s to %s\r",
					snc->chromosome_names[first_skipped],snc->chromosome_names[i]);
			continue;
		}
		printf("\nFinding Structural Variants in Chromosome %s\n",snc->chromosome_names[i]);
		first_skipped = i+1;
//		calculate_GC_histogram(in_bams, snc, i);

		short *depth_array = NULL;
		double global_molecule_mean = 0;
		double global_molecule_std = 0;
		regions[i] = recover_molecules(reads[i]->concordants);
		if( svs_to_find == SV_INVERSION || svs_to_find == SV_DUPLICATION || svs_to_find == SV_INVERTED_DUPLICATION){
			depth_array = make_molecule_depth_array(regions[i],snc,i);
			global_molecule_mean = make_global_molecule_mean(depth_array,snc,i);
			global_molecule_std = make_global_molecule_std_dev(depth_array,snc,i,global_molecule_mean);
			global_molecule_std = MIN(global_molecule_std,global_molecule_mean/2);
			printf("Global Molecule mean is: %lf\nGlobal Standard Deviation is: %lf\n",global_molecule_mean,global_molecule_std);
		}
		filter_molecules(regions[i]);

		qsort(regions[i]->items,regions[i]->size,sizeof(void*),interval_start_comp);  
		groups = group_overlapping_molecules(regions[i]);
		groups->rmv = &vector_free;

		CLONE_MEAN = molecule_mean(regions[i]);
		CLONE_STD_DEV = molecule_group_std(groups,CLONE_MEAN);

		vector_free(groups);

		qsort(regions[i]->items,regions[i]->size,sizeof(void*),barcode_comp);  

		printf("Chromosome %s\tClone Count %zu\tClone Mean %lf\tCLone Std %lf\n",snc->chromosome_names[i],regions[i]->size,CLONE_MEAN,CLONE_STD_DEV);

		vector_t *split_molecules = discover_split_molecules(regions[i]);
		vector_free(regions[i]);
		variations[i] = find_svs(split_molecules,svs_to_find);
		vector_free(split_molecules);
		printf("%zu candidate variations are made\n",variations[i]->size);
		update_sv_supports(variations[i],
				reads[i],
				svs_to_find);
		variations[i]->REMOVE_POLICY = REMP_LAZY;


		for (k=0;k<variations[i]->size;k++){
			sv_t *sv = vector_get(variations[i],k);
			sv_fprint(logFile,i,sv);
			if(sv->supports[0] < 1 || sv->supports[1] < 1){
				fprintf(logFile,"removed\n");
				vector_remove(variations[i],k);
			}

		}
		vector_defragment(variations[i]);
		if( svs_to_find == SV_INVERSION){
			for(k=0;k< variations[i]->size;k++){
				sv_t *inv = vector_get(variations[i],k);
				if(sonic_is_gap(snc,snc->chromosome_names[i],
							inv->AB.start1,inv->CD.end2)){
					sv_fprint(logFile,i,vector_get(variations[i],k));
					vector_remove(variations[i],k);
				}
			}
			vector_defragment(variations[i]);
		}
		if( svs_to_find == SV_DUPLICATION ||  svs_to_find == SV_INVERTED_DUPLICATION){
			for(k=0;k<variations[i]->size;k++){	
				sv_t *dup = vector_get(variations[i],k);
				int start = 0,end = 0,target_start = 0,target_end =0;
				if(dup->orientation == DUP_FORW_COPY){

					target_start = dup->AB.start2;
					target_end = dup->CD.end2;
					if(svs_to_find == SV_DUPLICATION){
						start = dup->AB.start1;
						end = dup->CD.end1;
					}
					else if(svs_to_find == SV_INVERTED_DUPLICATION){
						start = dup->CD.start1;
						end = dup->AB.end1;

					}
				}else if(dup->orientation == DUP_BACK_COPY){

					if(svs_to_find == SV_DUPLICATION){
						start = dup->AB.start2;
						end = dup->CD.end2;
					}else if(svs_to_find == SV_INVERTED_DUPLICATION){
						start = dup->CD.start2;
						end = dup->AB.end2;
					}
					target_start = dup->AB.start1;
					target_end = dup->CD.end1;
				}
				if(target_start > target_end){
					int temp = target_start;
					target_start = target_end;
					target_end = temp;
				}
				sonic_repeat *rep = sonic_is_mobile_element(snc,snc->chromosome_names[i],start,end,VALOR_MOBILE_ELEMENTS);
				if(rep!=NULL){
					int rep_start = MAX(rep->repeat_start,start);
					int rep_end = MIN(rep->repeat_end,end);
					if ( rep_end - rep_start < 0.6 * (end-start)){
						rep=NULL;
					}
				}
				int is_ref_dup_source = sonic_is_segmental_duplication(snc,snc->chromosome_names[i],start,end);
				int is_ref_dup_target = sonic_is_segmental_duplication(snc,snc->chromosome_names[i],target_start,target_end);

				int is_ref_gap_source = sonic_is_gap(snc,snc->chromosome_names[i],start,end);
				int is_ref_gap_target = sonic_is_gap(snc,snc->chromosome_names[i],target_start,target_end);
				int is_ref_sat_source = sonic_is_satellite(snc,snc->chromosome_names[i],start,end);
				int is_ref_sat_target = sonic_is_satellite(snc,snc->chromosome_names[i],target_start,target_end);	
				int does_cnv_support_dup = get_depth_region(depth_array,start,end) > global_molecule_mean + 1 * global_molecule_std;
				if( (is_ref_dup_source && is_ref_dup_target) || !does_cnv_support_dup   || is_ref_gap_source || is_ref_gap_target || is_ref_sat_source || is_ref_sat_target){
					fprintf(logFile,"%d %d %d %d %d %lf\n",does_cnv_support_dup, is_ref_dup_source, is_ref_dup_target, start, end, get_depth_region(depth_array,start,end));
					sv_fprint(logFile,i,vector_get(variations[i],k));
					vector_remove(variations[i],k);
				}
			}
				
			vector_defragment(variations[i]);
		}
		printf("%zu candidate variations are left after filtering\n",variations[i]->size);
#if FILTER1XK //If number of variations are more than MAX_INVERSIONS_IN_GRAPH, do random selection
		if(variations[i]->size > MAX_INVERSIONS_IN_GRAPH){
			srand(time(NULL));
			variations[i]->REMOVE_POLICY=REMP_LAZY;
			for(k=0;k<variations[i]->size;k++){
				if(rand()%variations[i]->size>MAX_INVERSIONS_IN_GRAPH){
					vector_remove(variations[i],k);
				}
			}
			vector_defragment(variations[i]);

			printf("%zu candidate variations are left after random selection to decrease size\n",variations[i]->size);
		}
#endif

		variations[i]->REMOVE_POLICY=REMP_SORTED;
		FILE *fii = fopen(OUT_DIR"/fbed.bedpe","w+");
		for(k=0;k<variations[i]->size;k++){
			sv_t *sv = vector_get(variations[i],k);
			sv_fprint(fii,i,sv);
		}
		fclose(fii);
		graph_t *sv_graph = make_sv_graph(variations[i]);
		printf("Finding Sv Clusters\n\n");
		clusters[i] = vector_init(sizeof(sv_cluster),50);
		vector_set_remove_function(clusters[i],&sv_cluster_destroy);



		graph_trim(sv_graph);

		vector_t *comps = g_dfs_components(sv_graph);

		if(comps->size == 0){continue;}

		for(k=0;k<comps->size;k++){
			vector_t *garbage = vector_get(comps,k);
			qsort(garbage->items, garbage->size, sizeof(sv_t *),&sv_comp);

/*
			if(comp->size < 16){break;}
			printf("\rCurrent Component Size %zu\t\t\t\t\n",comp->size);
			sv_cluster *svc = getMem(sizeof(sv_cluster));
			svc->items = vector_init(sizeof(sv_t),comp->size);
			VALOR_LOG("Cluster %d\n",k);
			sv_t *sv_to_add = (sv_t *)vector_get(comp,0);



			svc->supports[0] = sv_to_add->supports[0];
			svc->supports[1] = sv_to_add->supports[1];

			svc->break_points = sv_reduce_breakpoints(sv_to_add);	
			splitmolecule_t *scl_tmp;
			vector_t *garbage = vector_init(sizeof(sv_t),16);
			int ccc = 0;
			for(j=1;j<comp->size;j++){

				sv_to_add = (sv_t *)vector_get(comp,j);
				scl_tmp = sv_reduce_breakpoints(sv_to_add);
				if(!interval_pair_overlaps(svc->break_points,scl_tmp,CLONE_MEAN)){printf("\rCant Overlap %d items",++ccc);fflush(stdout);vector_put(garbage,sv_to_add);continue;}
				interval_pair_intersect(svc->break_points,scl_tmp);
				sv_fprint(logFile,i,sv_to_add);
				vector_put(svc->items, sv_to_add);
				svc->supports[0]+=sv_to_add->supports[0];
				svc->supports[1]+=sv_to_add->supports[1];
				splitmolecule_destroy(scl_tmp);
			}
			printf("\n");
			vector_put(clusters[i],svc);
//Now Explore Non Overlapping SV's
*/
			if(garbage->size < 4){continue;}
			graph_t *garbage_graph = sv_graph;//= make_sv_graph(garbage);
			printf("Component Size %zu\n",garbage->size);
			int iteration_no = 0;
			while(garbage_graph->number_of_items > 2){
				printf("graph size :%zu\n",garbage_graph->number_of_items);
				clique_t *c = clique_find_clique(garbage_graph,garbage,0,QCLIQUE_LAMBDA,QCLIQUE_GAMMA);

				if(c==NULL||c->v_prime<=0){clique_free(c);break;}
				sv_cluster *svc_garbage = sv_cluster_make(c);
				clique_free(c);
				
				if(svc_garbage==NULL){break;}
				if(svc_garbage->items->size < 4){
					sv_graph_reset(garbage_graph);
					sv_cluster_graph_fix(svc_garbage,garbage,sv_graph);
					sv_cluster_destroy(svc_garbage);
					continue;
				}
				printf("Found a clique of size %zu, Fixing Graph Now\n",svc_garbage->items->size);
				

				sv_graph_reset(garbage_graph);

				sv_cluster_graph_fix(svc_garbage,garbage,garbage_graph);



				vector_put(clusters[i],svc_garbage);
				iteration_no++;
			}
//			vector_free(garbage);
		//	graph_free(garbage_graph);
		}

		vector_free(comps);
/*

		while(sv_graph->number_of_items > 0){
			clique_t *c = clique_find_clique(sv_graph,0,QCLIQUE_LAMBDA,QCLIQUE_GAMMA);
			if(c==NULL||c->v_prime<=0 || c->items->size < 4){clique_free(c);break;}
			sv_cluster *svc = sv_cluster_make(c);
			clique_free(c);
			if(svc->items->size < 1){sv_cluster_graph_fix(svc,sv_graph);sv_cluster_destroy(svc);continue;}
			if(svc==NULL){break;}

			printf("Found a clique of size %zu, Fixing Graph Now\n",svc->items->size);
			sv_cluster_graph_fix(svc,sv_graph);

			//	graph_trim(sv_graph);

			vector_put(clusters[i],svc);
		}
*/
		qsort(clusters[i]->items, clusters[i]->size, sizeof( void*), cluster_comp);
		printf("Printing Clusters Now\n");
		for(j=0;j<clusters[i]->size;j++){
			sv_cluster *svc = vector_get(clusters[i],j);

			if(svc->items->size < 16){continue;}
			sv_t *first = vector_get(svc->items,0);


			double mean_depth = 0;
			if(svs_to_find != SV_INVERSION){
				if(first->orientation == DUP_FORW_COPY){
					mean_depth = get_depth_region(depth_array,svc->break_points->start1,svc->break_points->end1);
				}else if(first->orientation == DUP_BACK_COPY){
					mean_depth = get_depth_region(depth_array,svc->break_points->start1,svc->break_points->end1);
				}
				if(mean_depth < global_molecule_mean + 1.25 * global_molecule_std){ continue;}
			}else{
				mean_depth = get_depth_region(depth_array,first->AB.end1,first->CD.start1)/2 + get_depth_region(depth_array,first->AB.end2,first->CD.start2)/2;
			}

			fprintf(outbedfile,"%s\t%d\t%d\t%s\t%d\t%d\t%s\t%zu\t%d\t%lf\n",
					snc->chromosome_names[i],
					svc->break_points->start1,
					svc->break_points->end1,
					snc->chromosome_names[i],
					svc->break_points->start2,
					svc->break_points->end2,
					sv_type_name(svs_to_find),
					svc->items->size,
					svc->supports[0]+svc->supports[1],
					mean_depth       
			);
		}
		fflush(outbedfile);
#if DEVELOPMENT_
		for(j=0;j<clusters[i]->size;j++){
			sv_cluster *svc = vector_get(clusters[i],j);

			int is_dup =  sonic_is_segmental_duplication(snc,snc->chromosome_names[i],svc->break_points->start1,svc->break_points->end2);
			int is_gap =  sonic_is_gap(snc,snc->chromosome_names[i],svc->break_points->start1,svc->break_points->end2);

			fprintf(filecul,"Chromosome:%s\nCluster Size: %zu\nMM Support:%d\nPP Support:%d\nBreakpoints: %d\t%d\t%d\t%d\t%s\t%s\t\n",
					snc->chromosome_names[i], svc->items->size, svc->supports[1], svc->supports[0],svc->break_points->start1,
					svc->break_points->end1,svc->break_points->start2,svc->break_points->end2
					,is_dup?"On SegDup":"",is_gap?"On Gap":"");
			for(k=0;k<svc->items->size;k++){
				fprintf(filecul,"\t");
				sv_fprint(filecul,i,vector_get(svc->items,k));
			}
		}

		fclose(filecul);
#endif
		free(depth_array);
		destroy_bams(reads[i]);
		vector_free(variations[i]);
		vector_free(clusters[i]);
		graph_free(sv_graph);


		//		free(in_bams->read_depth[i]);
	}
	printf("\n");

	freeMem(in_bams, sizeof(bam_info));

	freeMem(clusters, sizeof(vector_t *) * snc->number_of_chromosomes);
	freeMem(variations,sizeof(vector_t *) * snc->number_of_chromosomes);
	free(regions);
	free(reads);
	free_sonic(snc);
	fclose(logFile);
	fclose(outbedfile);
	return 0;
}

