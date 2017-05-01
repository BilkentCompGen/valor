
#include "sonic_interval.h"
#include "sonic.h"


int bed_comp( const void* p1, const void* p2){
  sonic_bed_line a;
  sonic_bed_line b;
  int chrom_name_comp;
  
  a = * ( (sonic_bed_line *) p1);
  b = * ( (sonic_bed_line *) p2);

  chrom_name_comp = strcmp(a.chromosome, b.chromosome);

  if (chrom_name_comp != 0)
    return chrom_name_comp;

  return a.start - b.start;
 
}


int count_bed_chromosome_lines(FILE *bed_file, char *chromosome)
{
	int number_of_lines;
	char line[MAX_LENGTH];
	char *return_value;
	int return_value_int;
	char this_chromosome[MAX_LENGTH];
	
	number_of_lines = 0;
	while (!feof(bed_file)){
   	        return_value_int = fscanf(bed_file, "%s", this_chromosome);
		if (feof(bed_file) || return_value_int == 0)
			break;
		if (line[0] != 0){
 		        return_value = fgets(line, MAX_LENGTH, bed_file);
			if (return_value == NULL){
			  exit(EXIT_SONIC);
	  }

			if (!strcmp(this_chromosome, chromosome))
			        number_of_lines++;
		}
	}

	return number_of_lines;
}

/* sonic_interval *sonic_intersect(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end, sonic_interval_type interval_type, float intersection_fraction){ */
sonic_interval *sonic_intersect(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end, sonic_interval_type interval_type){


  int start;
  int end;
  int med;
  int interval_count;

  sonic_interval *this_interval_list;
  
  int chromosome_index;

  /*
  if (intersection_fraction < 0.0 || intersection_fraction > 1.0){
    fprintf(stderr, "[SONIC] Intersection fraction value should be in [0,1] interval. Fraction=0 assumes any intersection.\n");
    exit (EXIT_SONIC);
    }*/


  chromosome_index = sonic_find_chromosome_index(this_sonic->chromosome_names, this_chromosome, this_sonic->number_of_chromosomes);

  if (chromosome_index == -1)
    return NULL;
  
  switch (interval_type){
  case GAP:
    this_interval_list = this_sonic->gaps[chromosome_index];
    interval_count = this_sonic->number_of_gaps_in_chromosome[chromosome_index];
    break;
  case DUP:
    this_interval_list = this_sonic->dups[chromosome_index];
    interval_count = this_sonic->number_of_dups_in_chromosome[chromosome_index];
    break;
  case REP:
    this_interval_list = this_sonic->reps[chromosome_index];
    interval_count = this_sonic->number_of_repeats_in_chromosome[chromosome_index];
    break;
  default:
    return NULL;
  }
  
  start = 0;
  end = interval_count - 1;

  med = (start + end) / 2;

    
  while (1){

    //printf ("Search %d-%d-%d [%d]\n", start, med, end, interval_count);
    if (start > end)
      return NULL;

    if (sonic_this_interval_intersects(pos_start, pos_end, this_interval_list[med].start, this_interval_list[med].end))
      return &this_interval_list[med];

    /* no hit. search is exhausted */
    if (start == med || end == med){
      if (sonic_this_interval_intersects(pos_start, pos_end, this_interval_list[start].start, this_interval_list[start].end))
	return &this_interval_list[start];
      else if (sonic_this_interval_intersects(pos_start, pos_end, this_interval_list[end].start, this_interval_list[end].end))
	return &this_interval_list[end];      
      return NULL;
    }

    /*
    else if (start == med)
      med = end;
    else if (end == med)
      med = start;
    */
    
    /* no hit, search left half */
    else if (pos_start < this_interval_list[med].start){
      end = med;
      med = (start + end) / 2;
    }

    /* no hit, search right half */
    else {
      start = med;
      med = (start + end) / 2;      
    }
      
  }

  return NULL;
}


void sonic_print_interval(sonic_interval *this_interval){

  if (this_interval == NULL){
    fprintf(stdout, "[SONIC_INTERVAL] Not found.\n");
    return;
  }
  fprintf (stdout, "[SONIC_INTERVAL]\n");
  fprintf (stdout, "\tStart: %d\tEnd:%d\n", this_interval->start, this_interval->end);
  if (this_interval->repeat_item != NULL){
    fprintf (stdout, "\t[SONIC_INTERVAL_REPEAT]\n");
    fprintf (stdout, "\t\tType: %s\tClass: %s\n", this_interval->repeat_item->repeat_type, this_interval->repeat_item->repeat_class);
  }
  
}

int sonic_this_interval_intersects(int pos_start, int pos_end, int start, int end){
      /* all in */
    if (pos_start >= start && pos_end < end)
      return 1;

    /* all cover */
    else if (pos_start <= start && pos_end > end)
      return 1;

    /* left */
    else if (pos_start <= start && pos_end >= start)
      return 1;

    /* right */
    else if (pos_start <= end && pos_end > end)
      return 1;

    /* no hit.  */

    return 0;
}


int sonic_is_satellite(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end){
  sonic_interval *this_interval;
  char *is_satellite;

  /* potential problem here for general-case repeats within repeats */
  this_interval = sonic_intersect(this_sonic, this_chromosome, pos_start, pos_end, REP);

  if (this_interval == NULL)
    return 0;

  is_satellite = strstr(this_interval->repeat_item->repeat_class, "Satel");
  if (is_satellite != NULL)
    return 1;

  return 0;
    
  
}

int sonic_is_segmental_duplication(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end){

  sonic_interval *this_interval;

  this_interval = sonic_intersect(this_sonic, this_chromosome, pos_start, pos_end, DUP);

  if (this_interval == NULL)
    return 0;

  return 1;
    
}

int sonic_is_gap(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end){

  sonic_interval *this_interval;

  this_interval = sonic_intersect(this_sonic, this_chromosome, pos_start, pos_end, GAP);

  if (this_interval == NULL)
    return 0;

  return 1;

}

sonic_repeat *sonic_is_mobile_element(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end, char *mei_string){

  sonic_interval *this_interval;

  char *tok;
  char str[1024];
  
  this_interval = sonic_intersect(this_sonic, this_chromosome, pos_start, pos_end, REP);

  if (this_interval == NULL)
    return NULL;

  
  strcpy(str, mei_string);
  
  tok = strtok(str, ":");
  
  while (tok != NULL){
    if (strstr(this_interval->repeat_item->repeat_type, tok) != NULL){
      return this_interval->repeat_item;
    }
    
    tok = strtok(NULL, ":");
  }
    
  return NULL;
  
}


float sonic_get_gc_content(sonic *this_sonic, char *this_chromosome, int pos_start, int pos_end){

  
  int chromosome_index;
  int window_count;
  float gc_content;
  int start_gc;
  
  chromosome_index = sonic_find_chromosome_index(this_sonic->chromosome_names, this_chromosome, this_sonic->number_of_chromosomes);

  if (chromosome_index == -1)
    return 0.0;
  
  window_count = 0;
  gc_content = 0.0;

  start_gc = pos_start;

  while (start_gc < pos_end){
    printf("gc content window %s-%d. Index: %d. content %f\n", this_chromosome, start_gc, (start_gc/SONIC_GC_WINDOW), (float) this_sonic->chromosome_gc_profile[chromosome_index][start_gc / SONIC_GC_WINDOW]);
    window_count++;
    gc_content += (float) this_sonic->chromosome_gc_profile[chromosome_index][start_gc / SONIC_GC_WINDOW];
    start_gc += SONIC_GC_WINDOW;
  }

  if (gc_content != 0.0)
    return (gc_content / window_count);
  else
    return 0.0;
}
