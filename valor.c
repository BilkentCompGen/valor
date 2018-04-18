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
#include "vector.h"
#include "clique.h"
#include "cluster.h"
#include "cnv.h"
#include "sonic/sonic.h"
#include "graph.h"
#include "progress.h"
#include <omp.h>
int CUR_CHR = -1;
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
	vector_t **regions; // Vector of Inverval_10X
	vector_t *groups;   // Vector of Interval_10X
	vector_t **variations; //Vector of sv_t
	vector_t **clusters;   //Vector of cluster_t




	mkdir(params->outprefix, 0755 );
	char *out_file_path = malloc((strlen("/predicted_svs.bedpe")+strlen(params->outprefix)+1)*sizeof(char));
	sprintf(out_file_path,"%s/predicted_svs.bedpe",params->outprefix);
	FILE *outbedfile = fopen(out_file_path,"w+");
	
	time( &rawtime);

	timeinfo = localtime( &rawtime);


	printf("\n\nVALOR: Variation with LOng Range\n");
	printf("Version: %s\n", VALOR_VERSION);
	printf("Build Date: %s\n",BUILD_DATE);
	printf("Output Directory: %s\n",params->outprefix);
	printf("Logfile name: %s\n", params->logfile);


	sonic *snc = sonic_load(params->sonic_file);
	printf("Reading Bam file: %s\n", bamname);

	bam_info *in_bams = get_bam_info(snc);
	logFile = safe_fopen(params->logfile,"w+");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n",
			timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday); 
	in_bams->sample_name = NULL;


	bam_stats *stats = calculate_bam_statistics(in_bams, bamname, READ_SAMPLE_SIZE);

//	bam_vector_pack **reads = read_10X_bam(in_bams,bamname,snc);
	bam_vector_pack **reads = malloc(sizeof(bam_vector_pack) * snc->number_of_chromosomes);//read_10X_bam(in_bams,bamname,snc);

	regions = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	variations = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	clusters = getMem(sizeof(vector_t *) * snc->number_of_chromosomes);
	int first_skipped = 0;
	for( i = 0; i < 24 /*snc->number_of_chromosomes*/;i++){

		CUR_CHR = i;
		reads[i] = read_10X_chr(in_bams,bamname,snc,i,stats);
		if(reads[i]->concordants->size == 0){
			printf("No Reads for Chromosome %s %s %s\r",
				snc->chromosome_names[first_skipped],
				(i-first_skipped==1?"and":"to"),
				snc->chromosome_names[i]);
			continue;
		}
		printf("\nFinding Structural Variants in Chromosome %s\n",snc->chromosome_names[i]);
		first_skipped = i+1;

		printf("Recovering Split Molecules..\n");
		regions[i] = recover_molecules(reads[i]->concordants);

		in_bams->depths[i] = make_molecule_depth_array(regions[i],snc,i);
		
		in_bams->depth_mean[i] = make_global_molecule_mean(in_bams->depths[i],snc,i);
		in_bams->depth_std[i] = make_global_molecule_std_dev(in_bams->depths[i],snc,i,in_bams->depth_mean[i]);
		in_bams->depth_std[i] = MIN(in_bams->depth_std[i],in_bams->depth_mean[i]/2);
		printf("Global Molecule depth mean: %lf\nGlobal Molecule Depth Standard Deviation: %lf\n",
		in_bams->depth_mean[i],in_bams->depth_std[i]);


		filter_molecules(regions[i],snc,i);

		qsort(regions[i]->items,regions[i]->size,sizeof(void*),interval_start_comp);  
		groups = group_overlapping_molecules(regions[i]);
		groups->rmv = &vector_free;

		CLONE_MEAN = molecule_mean(regions[i]);
		CLONE_STD_DEV = molecule_group_std(groups,CLONE_MEAN);

		vector_free(groups);

		qsort(regions[i]->items,regions[i]->size,sizeof(void*),barcode_comp);  

		printf("Molecule Count: %zu\tMolecule Mean: %lf\tMolecule std-dev: %lf\n",regions[i]->size,CLONE_MEAN,CLONE_STD_DEV);

		vector_t *split_molecules = discover_split_molecules(regions[i]);
		vector_free(regions[i]);
		
		printf("Matching Split Molecules\n");
		variations[i] = find_svs(split_molecules,svs_to_find);
		vector_free(split_molecules);
		printf("%zu candidate variations are made\n",variations[i]->size);
		update_sv_supports(variations[i],
				reads[i],
				svs_to_find);
		variations[i]->REMOVE_POLICY = REMP_LAZY;

		vector_filter(variations[i],sv_is_proper);

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

		graph_t *sv_graph = make_sv_graph(variations[i]);
		printf("Finding Sv Clusters\n\n");
		clusters[i] = vector_init(sizeof(sv_cluster),50);
		vector_set_remove_function(clusters[i],&sv_cluster_destroy);



		graph_trim(sv_graph);

		vector_t *comps = g_dfs_components(sv_graph);

		if(comps->size == 0){continue;}
		
		for(k=0;k<comps->size;k++){
			vector_t *garbage = vector_get(comps,k);
			size_t initial_size = garbage->size;
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
			if(garbage->size < 16){continue;}
			graph_t *garbage_graph = sv_graph;//= make_sv_graph(garbage);
			int iteration_no = 0;
			while(garbage_graph->number_of_items > 2){
				if(garbage->size < initial_size /4){break;}

				clique_t *c = clique_find_clique(garbage_graph,garbage,0,QCLIQUE_LAMBDA,QCLIQUE_GAMMA);

				if(c==NULL||c->v_prime<=0){clique_free(c);break;}
				sv_cluster *svc_garbage = sv_cluster_make(c);
				clique_free(c);
				
				if(svc_garbage==NULL){break;}
				if(svc_garbage->items->size < 16){
					sv_graph_reset(garbage_graph);
					sv_cluster_graph_fix(svc_garbage,garbage,sv_graph);
					sv_cluster_destroy(svc_garbage);
					continue;
				}
	

				sv_graph_reset(garbage_graph);

				sv_cluster_graph_fix(svc_garbage,garbage,garbage_graph);



				vector_put(clusters[i],svc_garbage);
				iteration_no++;
			}
//			vector_free(garbage);
		//	graph_free(garbage_graph);
		}

		printf("Clustering is finished, found %zu variant clusters\n",clusters[i]->size);

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
		printf("Printing Variant calls\n");
		for(j=0;j<clusters[i]->size;j++){
			sv_cluster *svc = vector_get(clusters[i],j);

			if(svc->items->size < 16){continue;}
			sv_t *first = vector_get(svc->items,0);


			double mean_depth = 0;
			if(svs_to_find == SV_DUPLICATION || svs_to_find == SV_INVERTED_DUPLICATION){
				if(first->orientation == DUP_FORW_COPY){
					mean_depth = get_depth_region(in_bams->depths[i],svc->break_points->start1,svc->break_points->end1);
				}else if(first->orientation == DUP_BACK_COPY){
					mean_depth = get_depth_region(in_bams->depths[i],svc->break_points->start1,svc->break_points->end1);
				}
				if(mean_depth < in_bams->depth_mean[i] + 1.25 * in_bams->depth_std[i]){ continue;}
			}else{
				mean_depth = get_depth_region(in_bams->depths[i],first->AB.end1,first->CD.start1)/2 + get_depth_region(in_bams->depths[i],first->AB.end2,first->CD.start2)/2;
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
		destroy_bams(reads[i]);
		vector_free(variations[i]);
		vector_free(clusters[i]);
		graph_free(sv_graph);


		free(in_bams->depths[i]);
		//		free(in_bams->read_depth[i]);
	}
	printf("\n");

	for(i=0;i<snc->number_of_chromosomes;i++){

	}
	free(in_bams->depth_mean);
	free(in_bams->depth_std);
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

