/*
  SONIC - Some Organism's Nucleotide Information Container
*/

#include "sonic.h"
#include "sonic_interval.h"
#include "sonic_reference.h"

int sonic_build(char *ref_genome, char *gaps, char *reps, char *dups, char *info, char *sonic)
{
  FILE *ref_file;
  FILE *ref_index;
  FILE *gaps_file;
  FILE *reps_file;
  FILE *dups_file;
  gzFile sonic_file;
  int line_count;
  int sonic_magic = SONIC_MAGIC;
  int return_value;
  int chrom_name_length;
  char ref_genome_index[MAX_LENGTH];
  int i;
  time_t sonic_build_time;
  int info_length;
  
  int number_of_chromosomes;
  int *chromosome_lengths;
  char **chromosome_names;
  
  sonic_bed_line *bed_entry;

  sonic_build_time = time(NULL);
  sonic_mem_usage = 0;
  ref_file = sonic_fopen(ref_genome, "r");

  sprintf(ref_genome_index, "%s.fai", ref_genome);

  ref_index = sonic_fopen(ref_genome_index, "r");
  gaps_file = sonic_fopen(gaps, "r");
  reps_file = sonic_fopen(reps, "r");
  dups_file = sonic_fopen(dups, "r");

  sonic_file = sonic_fopen_gz(sonic, "w");

  return_value = gzwrite(sonic_file, &sonic_magic, sizeof(sonic_magic));

  if (return_value == 0){
    fprintf(stderr, "Cannot write magic number to SONIC file.\n");
    return EXIT_FILE_OPEN_ERROR;
  }

  
  return_value = gzwrite(sonic_file, &sonic_build_time, sizeof(sonic_build_time));

  if (info == NULL)
    sonic_set_str(&info, ref_genome);

  
  info_length = strlen(info);
  return_value = gzwrite(sonic_file, &info_length, sizeof(info_length));
  return_value = gzwrite(sonic_file, info, info_length);
  
  number_of_chromosomes = get_number_of_chromosomes(ref_index);
  rewind(ref_index);

  return_value = get_chromosome_info(ref_index, number_of_chromosomes, &chromosome_lengths, &chromosome_names);
  if (return_value == RETURN_ERROR){
    fprintf(stderr, "Cannot read the reference genome index properly. Check if %s is a valid index file.\n", ref_genome_index);
    return EXIT_SONIC;
  }

  fclose(ref_index);

  return_value = gzwrite(sonic_file, &number_of_chromosomes, sizeof(number_of_chromosomes));  /* write number of chromosomes */ 
  fprintf(stderr, "Number of chromosomes: %d\n", number_of_chromosomes);
  
  for (i=0; i < number_of_chromosomes; i++){  /* write chromosome names and lengths */ 
    /* fprintf(stderr, "Chromosome name: %s, length: %d\n", chromosome_names[i], chromosome_lengths[i]); */
    chrom_name_length = strlen(chromosome_names[i]);
    return_value = gzwrite(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzwrite(sonic_file, chromosome_names[i], chrom_name_length);
    return_value = gzwrite(sonic_file, &chromosome_lengths[i], sizeof(int));      
  }


    
  fprintf( stderr, "Adding gap intervals to SONIC.\n");
  line_count = count_bed_lines(gaps_file);
  rewind(gaps_file);

  bed_entry = sonic_read_bed_file(gaps_file, line_count, 0);
  sonic_write_bed_entries(sonic_file, bed_entry, line_count, number_of_chromosomes, chromosome_names);  
  fclose(gaps_file);
  free(bed_entry);


  fprintf( stderr, "Adding segmental duplication intervals to SONIC.\n");
  line_count = count_bed_lines(dups_file);
  rewind(dups_file);

  bed_entry = sonic_read_bed_file(dups_file, line_count, 0);
  sonic_write_bed_entries(sonic_file, bed_entry, line_count, number_of_chromosomes, chromosome_names);
  fclose(dups_file);
  free(bed_entry);

  
  fprintf( stderr, "Adding repeats to SONIC.\n");
  line_count = count_bed_lines(reps_file);
  rewind(reps_file);

  bed_entry = sonic_read_bed_file(reps_file, line_count, 1);
  sonic_write_bed_entries(sonic_file, bed_entry, line_count, number_of_chromosomes, chromosome_names);
  fclose(reps_file);
  
  for (i=0; i < line_count; i++){
    if (bed_entry[i].repeat_item != NULL){
      free(bed_entry[i].repeat_item->repeat_type);
      free(bed_entry[i].repeat_item->repeat_class);
      free(bed_entry[i].repeat_item);
    }
  }
  
  free(bed_entry);	 


  /* gc profile here */

  sonic_write_gc_profile(sonic_file, ref_file, number_of_chromosomes, chromosome_names);
  
  for (i=0; i < number_of_chromosomes; i++)  /* free memory */
    free(chromosome_names[i]);
  free(chromosome_names);
  free(chromosome_lengths);
  

  gzclose(sonic_file);

  fprintf( stderr, "SONIC file %s is ready.\n", sonic);
  fprintf( stdout, "Memory usage: %0.2f MB.\n", sonic_get_mem_usage());
  return RETURN_SUCCESS;
  
}



sonic *sonic_load(char *sonic_file_name){

  gzFile sonic_file;
  int sonic_magic;
  int return_value;
  char chromosome[MAX_LENGTH];
  int chrom_name_length;
  int start, end;
  int number_of_chromosomes;
  int number_of_entries;
  int chromosome_length;
  int i, j;

  /* 
     RepeatMasker out file entries
  */


  char strand;
  char repeat_type[MAX_LENGTH];
  char repeat_class[MAX_LENGTH];
  int repeat_type_length;
  int repeat_class_length;
  int repeat_start, repeat_end;
  sonic *this_sonic;

  time_t sonic_build_time;
  struct tm *sonic_ptm;
  int info_length;
  char *info;
  
  sonic_mem_usage = 0;
  
  fprintf (stderr, "Loading SONIC file..\n");
		     
  sonic_file = sonic_fopen_gz(sonic_file_name, "r");

  
  return_value = gzread(sonic_file, &sonic_magic, sizeof(sonic_magic));

  if (return_value == 0){
    fprintf(stderr, "Cannot read the SONIC file.\n");
    return NULL;
  }

  if (sonic_magic != SONIC_MAGIC){
    fprintf(stderr, "Invalid SONIC file (%s). Please use a correctly created SONIC file.\n", sonic_file_name);
    return NULL;
  }


  return_value = gzread(sonic_file, &sonic_build_time, sizeof(sonic_build_time));

  sonic_ptm = gmtime(&sonic_build_time);
  
  return_value = gzread(sonic_file, &info_length, sizeof(info_length));
  info = (char *) sonic_get_mem(sizeof(char *) * (info_length+1));
  
  return_value = gzread(sonic_file, info, info_length);
  info[info_length] = 0;

  fprintf(stderr, "\nSONIC Info: %s\nBuilt in %s\n\n", info, asctime(sonic_ptm));
  
  sonic_free_mem(info, (info_length+1));
  
  return_value = gzread(sonic_file, &number_of_chromosomes, sizeof(number_of_chromosomes));  /* read number of chromosomes */ 
  fprintf(stderr, "Number of chromosomes: %d\n", number_of_chromosomes);

  this_sonic = alloc_sonic(number_of_chromosomes);
  
  for (i=0; i < number_of_chromosomes; i++){  /* read chromosome names and lengths */ 
    return_value = gzread(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzread(sonic_file, chromosome, chrom_name_length);
    chromosome[chrom_name_length] = 0;
    return_value = gzread(sonic_file, &chromosome_length, sizeof(int));
    sonic_set_str(&(this_sonic->chromosome_names[i]), chromosome);
    this_sonic->chromosome_gc_profile[i] = (char *) sonic_get_mem(sizeof(char ) * (chromosome_length / (SONIC_GC_WINDOW)));
    this_sonic->chromosome_lengths[i] = chromosome_length;
    this_sonic->genome_length += chromosome_length;
    /* fprintf(stderr, "Chromosome name: %s, length: %d\n", chromosome, chromosome_length); */
  }

  
  fprintf( stderr, "Loading gap intervals.\n");


  for (i=0; i < number_of_chromosomes; i++){
    return_value = gzread(sonic_file, &number_of_entries, sizeof(number_of_entries)); /* number of gaps in this chromosome */
    
    this_sonic->number_of_gaps_in_chromosome[i] = number_of_entries;
    if (number_of_entries != 0)
      this_sonic->gaps[i] = alloc_sonic_interval(number_of_entries, 0);
    else
      this_sonic->gaps[i] = NULL;
    
    for (j = 0; j < number_of_entries; j++){
      return_value = gzread(sonic_file, &start, sizeof(start));
      return_value = gzread(sonic_file, &end, sizeof(end));
      this_sonic->gaps[i][j].start = start;
      this_sonic->gaps[i][j].end = end;
      /* fprintf(stderr, "[gaps]\t%d\t%d\n", start, end); */
    }
  }    

  fprintf( stderr, "Loading duplication intervals.\n");

  for (i=0; i < number_of_chromosomes; i++){
    return_value = gzread(sonic_file, &number_of_entries, sizeof(number_of_entries)); /* number of gaps in this chromosome */
    /* fprintf(stderr, "Chromosome %d dups %d\n", i, number_of_entries); */
    
    this_sonic->number_of_dups_in_chromosome[i] = number_of_entries;
    if (number_of_entries != 0)  
      this_sonic->dups[i] = alloc_sonic_interval(number_of_entries, 0);
    else
      this_sonic->dups[i] = NULL;
    
    for (j = 0; j < number_of_entries; j++){
      return_value = gzread(sonic_file, &start, sizeof(start));
      return_value = gzread(sonic_file, &end, sizeof(end));
      this_sonic->dups[i][j].start = start;
      this_sonic->dups[i][j].end = end;
      /*fprintf(stderr, "[dups]\t%d\t%d\n", start, end); */
    }

  }    


  fprintf( stderr, "Loading repeats.\n");

  for (i=0; i < number_of_chromosomes; i++){
    return_value = gzread(sonic_file, &number_of_entries, sizeof(number_of_entries)); /* number of gaps in this chromosome */
    /*    fprintf(stderr, "Chromosome %d reps %d\n", i, number_of_entries); */

    this_sonic->number_of_repeats_in_chromosome[i] = number_of_entries;
    if (number_of_entries != 0)
      this_sonic->reps[i] = alloc_sonic_interval(number_of_entries, 1);
    else
      this_sonic->reps[i] = NULL;

    for (j = 0; j < number_of_entries; j++){
      return_value = gzread(sonic_file, &start, sizeof(start));
      return_value = gzread(sonic_file, &end, sizeof(end));

      return_value = gzread(sonic_file, &strand, sizeof(strand));

      return_value = gzread(sonic_file, &repeat_type_length, sizeof(repeat_type_length));
      return_value = gzread(sonic_file, repeat_type, repeat_type_length);
      return_value = gzread(sonic_file, &repeat_class_length, sizeof(repeat_class_length));
      return_value = gzread(sonic_file, repeat_class, repeat_class_length);
      repeat_type[repeat_type_length] = 0;
      repeat_class[repeat_class_length] = 0;
      return_value = gzread(sonic_file, &repeat_start, sizeof(repeat_start));
      return_value = gzread(sonic_file, &repeat_end, sizeof(repeat_end));
           	    
      this_sonic->reps[i][j].start = start;
      this_sonic->reps[i][j].end = end;

      this_sonic->reps[i][j].repeat_item->strand = strand;
      this_sonic->reps[i][j].repeat_item->repeat_start = repeat_start;
      this_sonic->reps[i][j].repeat_item->repeat_end = repeat_end;
      sonic_set_str(&(this_sonic->reps[i][j].repeat_item->repeat_class), repeat_class);
      sonic_set_str(&(this_sonic->reps[i][j].repeat_item->repeat_type), repeat_type);
      
      /* fprintf(stderr, "[reps]\t%d\t%d\t%c\t%s\t%s\t%d\t%d\n", start, end, strand, repeat_type, repeat_class, repeat_start, repeat_end); */

    }

  }    

  /* read GC profiles */

  
  sonic_read_gc_profile(sonic_file, this_sonic);
  
  gzclose(sonic_file);

  fprintf( stdout, "SONIC file loaded. Memory usage: %0.2f MB.\n", sonic_get_mem_usage());
  return this_sonic;

  
}

FILE* sonic_fopen( char* path, const char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;

	file = fopen( path, mode);  
	if( !file)
	{
		fprintf( stderr, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.\n", path, mode[0]=='w' ? "write" : "read");
		exit(EXIT_FILE_OPEN_ERROR);

	}
	return file;
}

gzFile sonic_fopen_gz( char* path, const char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	gzFile file;

	file = gzopen( path, mode);  
	if( !file)
	{
	        fprintf( stderr, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.\n", path, mode[0]=='w' ? "write" : "read");
		exit(EXIT_FILE_OPEN_ERROR);		
	}
	return file;
}

int count_bed_lines(FILE *bed_file)
{
	int number_of_lines;
	char line[MAX_LENGTH];
	char *return_value;

	number_of_lines = 0;
	while (!feof(bed_file)){
		return_value = fgets(line, MAX_LENGTH, bed_file);
		if (feof(bed_file) || return_value == NULL)
			break;
		if (line[0] != 0)
			number_of_lines++;
	}

	return number_of_lines;
}


int count_bed_chromosome_entries(sonic_bed_line *bed_entry, int number_of_entries, char *chromosome)
{
        int count;
	int i;
	int found_chromosome;

	/* assume sorted bed. This is correct since we ran qsort before using this */
	
	count = 0;
	found_chromosome = 0;
	
	for (i=0; i < number_of_entries; i++){
	  if (!strcmp(bed_entry[i].chromosome, chromosome)){
	    found_chromosome = 1;
	    count++;
	  }
	  else if (found_chromosome){
	    /* the chromosome was already found before, now something else comes up -- end of this chromosome. */
	    return count;
	  }
	}

	/* end of lines */
	return count;
}

void sonic_set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	if (source != NULL)
	{
		( *target) = ( char*) sonic_get_mem( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
		( *target) = NULL;
	}
}

void* sonic_get_mem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "[SONIC] Cannot allocate memory. Requested memory = %0.2f MB.\n", ( float) ( size / 1048576.0));
		exit( 0);
	}

	sonic_mem_usage += size;
	return ret;
}

void sonic_free_mem( void *ptr, size_t size)
{

	if ( ptr != NULL){
	  
	  free( ptr);
	  sonic_mem_usage -= size;
	  
	}
}

void sonic_write_bed_entries(gzFile sonic_file, sonic_bed_line *bed_entry, int line_count, int number_of_chromosomes, char **chromosome_names)
{
  int i, j;
  int return_value;
  int chromosome_found;
  int number_of_entries;
  int wrote;
  
  wrote = 0;

  
  for (i=0; i < number_of_chromosomes; i++){
    number_of_entries = count_bed_chromosome_entries(bed_entry, line_count, chromosome_names[i]);
    
    return_value = gzwrite(sonic_file, &number_of_entries, sizeof(number_of_entries)); /* number of gaps in this chromosome */
    chromosome_found = 0;
    for (j = 0; j < line_count; j++){
      if (!strcmp(bed_entry[j].chromosome, chromosome_names[i])){
	  chromosome_found = 1;
	  return_value = gzwrite(sonic_file, &bed_entry[j].start, sizeof(bed_entry[j].start));
	  return_value = gzwrite(sonic_file, &bed_entry[j].end, sizeof(bed_entry[j].end));
	  if (return_value == 0){
	    exit(EXIT_SONIC);
	  }
	  if (bed_entry[j].repeat_item != NULL)
	    sonic_write_repeat_item(sonic_file, bed_entry[j].repeat_item);
	  
	  wrote++;
	}
      else if (chromosome_found)
        break; /* this chromosome is finished */
    }    
  }

  fprintf(stderr, "Wrote %d entries.\n", wrote);
}

sonic_bed_line *sonic_read_bed_file(FILE *bed_file, int line_count, int is_repeat)
{

  int i;
  char chromosome[255];
  int start, end;
  sonic_bed_line *bed_entry;
  int return_value;
  char *return_value_char; 

  char skip_line[MAX_LENGTH];
  char sw_score[MAX_LENGTH];
  float perc_div, perc_del, perc_ins;
  char chrom_left[MAX_LENGTH];
  char strand[2];
  char repeat_type[MAX_LENGTH];
  char repeat_class[MAX_LENGTH];
  int repeat_start, repeat_end;
  char repeat_start_string[MAX_LENGTH];
  char repeat_end_string[MAX_LENGTH];
  char repeat_left[MAX_LENGTH];
  int repeat_id;
  

  bed_entry = (sonic_bed_line *) sonic_get_mem(sizeof(sonic_bed_line) * line_count);

  i = 0;

  if (!is_repeat){
    while (fscanf(bed_file, "%s\t%d\t%d\n", chromosome, &start, &end) > 0){
      strncpy(bed_entry[i].chromosome, chromosome, strlen(chromosome)+1);
      bed_entry[i].start = start;
      bed_entry[i].end = end;
      bed_entry[i].repeat_item = NULL;
      i++;
    }
  }

  else {

    /* assuming vanilla RepeatMasker .out files here -- concatenated */
    while (fscanf(bed_file, "%s", sw_score) > 0){
      
      if (!strcmp(sw_score, "SW")){
	return_value_char = fgets(skip_line, MAX_LENGTH, bed_file); /* skip the first header line */
	return_value_char = fgets(skip_line, MAX_LENGTH, bed_file); /* skip the second header line */
	return_value_char = fgets(skip_line, MAX_LENGTH, bed_file); /* skip the empty line */
	if (return_value_char == NULL){
	  exit(EXIT_SONIC);
	}

	continue;
      }
      
      /* data starts here */
      
      return_value = fscanf(bed_file, "%f%f%f%s%d%d%s%s%s%s%s%s%s%d\n", &perc_div, &perc_del, &perc_ins, chromosome,
			    &start, &end, chrom_left, strand, repeat_type, repeat_class, repeat_start_string, repeat_end_string, repeat_left, &repeat_id);
      
      if (return_value == 0){
	exit(EXIT_SONIC);
      }
      
      /* keep the necessary ones */

      strncpy(bed_entry[i].chromosome, chromosome, strlen(chromosome)+1);
      bed_entry[i].start = start;
      bed_entry[i].end = end;     
      bed_entry[i].repeat_item = (sonic_repeat *) sonic_get_mem(sizeof(sonic_repeat));

      if (strand[0] == '+'){
	bed_entry[i].repeat_item->strand = SONIC_STRAND_FWD;
	repeat_start = atoi(repeat_start_string);
	repeat_end   = atoi(repeat_end_string);
      }
      else{
      	bed_entry[i].repeat_item->strand = SONIC_STRAND_REV;
	repeat_start = atoi(repeat_left);
	repeat_end   = atoi(repeat_end_string);
      }

      bed_entry[i].repeat_item->repeat_type = NULL;
      bed_entry[i].repeat_item->repeat_class = NULL;
      
      sonic_set_str(&(bed_entry[i].repeat_item->repeat_type), repeat_type);
      sonic_set_str(&(bed_entry[i].repeat_item->repeat_class), repeat_class);

      bed_entry[i].repeat_item->repeat_start= repeat_start;
      bed_entry[i].repeat_item->repeat_end= repeat_end;
      
      i++;
    }
    
  }
  
  qsort(bed_entry, line_count, sizeof(sonic_bed_line), bed_comp); /* sort the bed entries */
  return bed_entry;
}


void sonic_write_repeat_item(gzFile sonic_file, sonic_repeat *repeat_item)
{
  int return_value;
  int repeat_type_length;
  int repeat_class_length;
  
  return_value = gzwrite(sonic_file, &(repeat_item->strand), sizeof(repeat_item->strand));
  repeat_type_length = strlen(repeat_item->repeat_type);
  repeat_class_length = strlen(repeat_item->repeat_class);
  return_value = gzwrite(sonic_file, &repeat_type_length, sizeof(repeat_type_length));
  return_value = gzwrite(sonic_file, repeat_item->repeat_type, repeat_type_length);
  return_value = gzwrite(sonic_file, &repeat_class_length, sizeof(repeat_class_length));
  return_value = gzwrite(sonic_file, repeat_item->repeat_class, repeat_class_length);

  return_value = gzwrite(sonic_file, &repeat_item->repeat_start, sizeof(repeat_item->repeat_start));
  return_value = gzwrite(sonic_file, &repeat_item->repeat_end, sizeof(repeat_item->repeat_end));
  
  if (return_value == 0){
    exit(EXIT_SONIC);
  }

}


sonic *alloc_sonic(int number_of_chromosomes)
{
  sonic *new_sonic;
  int i;
  
  new_sonic = (sonic *) sonic_get_mem (sizeof(sonic));


  new_sonic->chromosome_lengths = (int *) sonic_get_mem(sizeof(int) * number_of_chromosomes);
  new_sonic->number_of_gaps_in_chromosome = (int *) sonic_get_mem(sizeof(int) * number_of_chromosomes);
  new_sonic->number_of_dups_in_chromosome = (int *) sonic_get_mem(sizeof(int) * number_of_chromosomes);
  new_sonic->number_of_repeats_in_chromosome = (int *) sonic_get_mem(sizeof(int) * number_of_chromosomes);

  new_sonic->chromosome_names = (char **) sonic_get_mem(sizeof(char *) * number_of_chromosomes);
  new_sonic->chromosome_gc_profile = (char **) sonic_get_mem(sizeof(char *) * number_of_chromosomes);
  
  new_sonic->gaps = (sonic_interval **) sonic_get_mem(sizeof(sonic_interval *) * number_of_chromosomes);
  new_sonic->dups = (sonic_interval **) sonic_get_mem(sizeof(sonic_interval *) * number_of_chromosomes);
  new_sonic->reps = (sonic_interval **) sonic_get_mem(sizeof(sonic_interval *) * number_of_chromosomes);

  for (i = 0; i < number_of_chromosomes; i++){
    new_sonic->chromosome_names[i] = NULL;
  }
  
  new_sonic->number_of_chromosomes = number_of_chromosomes;
  new_sonic->genome_length = 0;

  new_sonic->last_chromosome_index = -1;
  return new_sonic;
}

sonic_interval *alloc_sonic_interval(int number_of_entries, int is_repeat)
{
  sonic_interval *new_sonic_interval;
  int i;
  
  new_sonic_interval = (sonic_interval *) sonic_get_mem (sizeof(sonic_interval) * number_of_entries);

  if (!is_repeat)
    new_sonic_interval->repeat_item = NULL;
  else{
    for (i = 0; i < number_of_entries; i++){
      new_sonic_interval[i].repeat_item = (sonic_repeat *) sonic_get_mem (sizeof(sonic_repeat));
      new_sonic_interval[i].repeat_item->repeat_type = NULL;
      new_sonic_interval[i].repeat_item->repeat_class = NULL;
      new_sonic_interval[i].repeat_item->repeat_start = 0;
      new_sonic_interval[i].repeat_item->repeat_end = 0;
      new_sonic_interval[i].repeat_item->strand = 0;
    }
  }

  return new_sonic_interval;
}

void free_sonic(sonic *this_sonic)
{

  int i;

  int number_of_chromosomes;

  number_of_chromosomes = this_sonic->number_of_chromosomes;

  for (i = 0; i < number_of_chromosomes; i++){
    sonic_free_mem(this_sonic->chromosome_names[i], strlen(this_sonic->chromosome_names[i]));
    sonic_free_mem(this_sonic->chromosome_gc_profile[i], sizeof(char *) * (this_sonic->chromosome_lengths[i] / (SONIC_GC_WINDOW)));
    free_sonic_interval(this_sonic->gaps[i], this_sonic->number_of_gaps_in_chromosome[i], 0);
    free_sonic_interval(this_sonic->dups[i], this_sonic->number_of_dups_in_chromosome[i], 0);
    free_sonic_interval(this_sonic->reps[i], this_sonic->number_of_repeats_in_chromosome[i], 1);
  }
  

  sonic_free_mem(this_sonic->chromosome_lengths, sizeof(int) * number_of_chromosomes);
  
  sonic_free_mem(this_sonic->number_of_gaps_in_chromosome, sizeof(int) * number_of_chromosomes);
  sonic_free_mem(this_sonic->number_of_dups_in_chromosome, sizeof(int) * number_of_chromosomes);
  sonic_free_mem(this_sonic->number_of_repeats_in_chromosome, sizeof(int) * number_of_chromosomes);

  
  sonic_free_mem(this_sonic->chromosome_names, sizeof(char *) * number_of_chromosomes);
  sonic_free_mem(this_sonic->chromosome_gc_profile, sizeof(char *) * number_of_chromosomes); 
  
  sonic_free_mem(this_sonic->gaps, sizeof(sonic_interval *) * number_of_chromosomes);
  sonic_free_mem(this_sonic->dups, sizeof(sonic_interval *) * number_of_chromosomes);
  sonic_free_mem(this_sonic->reps, sizeof(sonic_interval *) * number_of_chromosomes);
  
  sonic_free_mem(this_sonic, sizeof(sonic));
  
}

void free_sonic_interval(sonic_interval *this_sonic_interval, int number_of_entries, int is_repeat)
{

  int i;
  

  if (is_repeat){    
    for (i = 0; i < number_of_entries; i++){
      sonic_free_mem (this_sonic_interval[i].repeat_item->repeat_type, strlen(this_sonic_interval[i].repeat_item->repeat_type));
      sonic_free_mem (this_sonic_interval[i].repeat_item->repeat_class, strlen(this_sonic_interval[i].repeat_item->repeat_class));
      sonic_free_mem (this_sonic_interval[i].repeat_item, sizeof(sonic_repeat));
    }
  }

  sonic_free_mem (this_sonic_interval, sizeof(sonic_interval) * number_of_entries);

}


double sonic_get_mem_usage()
{
  return sonic_mem_usage/1048576.0;
}


