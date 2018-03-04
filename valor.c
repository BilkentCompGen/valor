#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/ioctl.h>  

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

#include "progress.h"
int CUR_CHR = -1;
FILE *logFile = NULL;
double CLONE_MEAN;
double CLONE_STD_DEV;

int main( int argc, char **argv){

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
	FILE *outbedfile = fopen(OUT_DIR"/predicted_inversions.bed","w+");

	time( &rawtime);
	timeinfo = localtime( &rawtime);

	logFile = safe_fopen("valor.log","w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n",
		 timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday); 

	in_bams = (bam_info *) getMem( sizeof(bam_info));
	in_bams->sample_name = NULL;
	bam_vector_pack *reads = read_10X_bam(in_bams,bamname,snc);

	regions = getMem(sizeof(vector_t *) * reads->count);
	variations = getMem(sizeof(vector_t *) * reads->count);
	clusters = getMem(sizeof(vector_t *) * reads->count);
	int first_skipped = 0;
	for( i=0; i<reads->count;i++){
#if DEVELOPMENT_
		char fileculname[255];
		sprintf(fileculname,OUT_DIR"clusters%d.out",i+1);
		FILE *filecul = fopen( fileculname,"w+");
		if(filecul==NULL){ fprintf(stderr,"Could not open %s\n",fileculname); return -1;}
#endif
		CUR_CHR = i;
		if(reads->concordants[i]->size == 0){
			if( i == first_skipped){
			printf("No Reads for Chromosome %s\r",snc->chromosome_names[first_skipped]);
			}else{
			printf("No Reads for Chromosome %s to %s\r",
				snc->chromosome_names[first_skipped],snc->chromosome_names[i]);
			}
			free(in_bams->read_depth[i]);
			continue;
		}
		printf("\nFinding Inversions in Chromosome %s\n",snc->chromosome_names[i]);
		first_skipped = i+1;
		calculate_GC_histogram(in_bams, snc, i);

		regions[i] = recover_molecules(reads->concordants[i]);
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
//		vector_free(regions[i]);
		variations[i] = find_svs(split_molecules,SV_INVERSION);
//		vector_free(split_molecules);
		printf("%zu candidate variations are made\n",variations[i]->size);
		update_sv_supports(variations[i],
			reads->mm_discordants[i],
			reads->pp_discordants[i],
			SV_INVERSION);
		variations[i]->REMOVE_POLICY = REMP_LAZY;
	
		for (k=0;k<variations[i]->size;k++){
			sv_t *sv = vector_get(variations[i],k);
			if(sv->supports[0] < 1 || sv->supports[1] < 1){
				vector_remove(variations[i],k);
			}

		}
		vector_defragment(variations[i]);
		for(k=0;k< variations[i]->size;k++){
			sv_t *inv = vector_get(variations[i],k);
			if(sonic_is_gap(snc,snc->chromosome_names[i],
				inv->AB.start1,inv->CD.end2)){
				vector_remove(variations[i],k);
			}
		}
		vector_defragment(variations[i]);

		printf("%zu candidate variations are left after read pair support filtering\n",variations[i]->size);
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
		for(k=0;k<variations[i]->size;k++){
			sv_fprint(logFile,i,vector_get(variations[i],k));
		}
		variations[i]->REMOVE_POLICY = REMP_SORTED;
		
		graph_t *sv_graph = make_sv_graph(variations[i]);

		sv_graph->hf = &sv_hf;
		graph_set_rem_function(sv_graph,&sv_destroy);
		clusters[i] = vector_init(sizeof(sv_cluster),50);
		vector_set_remove_function(clusters[i],&sv_cluster_destroy);
		printf("Graph Node #:%zu\n", sv_graph->number_of_items);
		
		while(sv_graph->number_of_items > 0){
			clique_t *c = clique_find_clique(sv_graph,0,QCLIQUE_LAMBDA,QCLIQUE_GAMMA);
			if(c==NULL||c->v_prime<=0){clique_free(c);break;}
			sv_cluster *svc = sv_cluster_make(c);
			clique_free(c);
			if(svc->items->size < 1){sv_cluster_graph_fix(svc,sv_graph);sv_cluster_destroy(svc);continue;}
			if(svc==NULL){break;}
			sv_cluster_graph_fix(svc,sv_graph);
		//	graph_trim(sv_graph);

			vector_put(clusters[i],svc);
		}
		qsort(clusters[i]->items, clusters[i]->size, sizeof( void*), cluster_comp);


		for(j=0;j<clusters[i]->size;j++){
			sv_cluster *svc = vector_get(clusters[i],j);
			fprintf(outbedfile,"%s\t%d\t%d\t%d\t%d\t%zu\t%d\n",
			snc->chromosome_names[i],
			svc->break_points->start1,
			svc->break_points->end1,
			svc->break_points->start2,
			svc->break_points->end2,
			svc->items->size,
			svc->supports[0]+svc->supports[1]
		
			);
		}
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

//		vector_free(variations[i]);
//		vector_free(clusters[i]);
//		graph_free(sv_graph);
		fclose(filecul);
#endif
//		free(in_bams->read_depth[i]);
	}
	printf("\n");

	freeMem(in_bams, sizeof(bam_info));

	freeMem(clusters, sizeof(vector_t *) * reads->count);
	freeMem(variations,sizeof(vector_t *) * reads->count);
	free(regions);
	free_sonic(snc);
	

	destroy_bams(reads);
	fclose(logFile);
	fclose(outbedfile);
	return 0;
}

