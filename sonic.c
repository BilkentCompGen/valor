/*
  SONIC - Some Organism's Nucleotide Information Container
*/

#include "sonic.h"


int make_sonic(char *ref_genome, char *gaps, char *reps, char *dups, char *sonic)
{
  FILE *ref_file;
  FILE *gaps_file;
  FILE *reps_file;
  FILE *dups_file;
  gzFile sonic_file;
  int line_count;
  int sonic_magic = SONIC_MAGIC;
  int return_value;
  char *return_value_char;
  char chromosome[MAX_LENGTH];
  int chrom_name_length;
  int start, end;


  /* 
     RepeatMasker out file entries
  */

  char skip_line[MAX_LENGTH];
  char sw_score[MAX_LENGTH];
  float perc_div, perc_del, perc_ins;
  char chrom_left[MAX_LENGTH];
  char strand[2];
  char repeat_type[MAX_LENGTH];
  char repeat_class[MAX_LENGTH];
  int repeat_type_length;
  int repeat_class_length;
  int repeat_start, repeat_end;
  char repeat_start_string[MAX_LENGTH];
  char repeat_end_string[MAX_LENGTH];
  char repeat_left[MAX_LENGTH];
  int repeat_id;
  
  
  ref_file = sonic_fopen(ref_genome, "r");
  gaps_file = sonic_fopen(gaps, "r");
  reps_file = sonic_fopen(reps, "r");
  dups_file = sonic_fopen(dups, "r");

  sonic_file = sonic_fopen_gz(sonic, "w");

  return_value = gzwrite(sonic_file, &sonic_magic, sizeof(sonic_magic));

  if (return_value == 0){
    fprintf(stderr, "Cannot write magic number to SONIC file.\n");
    return EXIT_FILE_OPEN_ERROR;
  }

  line_count = count_bed_lines(gaps_file);
  rewind(gaps_file);

  fprintf( stderr, "Adding %d gap intervals to SONIC.\n", line_count);
  return_value = gzwrite(sonic_file, &line_count, sizeof(line_count));
  
  while (fscanf(gaps_file, "%s\t%d\t%d\n", chromosome, &start, &end) > 0){
    chrom_name_length = strlen(chromosome);
    return_value = gzwrite(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzwrite(sonic_file, chromosome, chrom_name_length);
    return_value = gzwrite(sonic_file, &start, sizeof(start));
    return_value = gzwrite(sonic_file, &end, sizeof(end));
  }

  fclose(gaps_file);
 
  line_count = count_bed_lines(dups_file);
  rewind(dups_file);

  return_value = gzwrite(sonic_file, &line_count, sizeof(line_count));

  fprintf( stderr, "Adding %d dup intervals to SONIC.\n", line_count);
  while (fscanf(dups_file, "%s\t%d\t%d\n", chromosome, &start, &end) > 0){
    chrom_name_length = strlen(chromosome);
    return_value = gzwrite(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzwrite(sonic_file, chromosome, chrom_name_length);
    return_value = gzwrite(sonic_file, &start, sizeof(start));
    return_value = gzwrite(sonic_file, &end, sizeof(end));
  }

  fclose(dups_file);

  fprintf( stderr, "Adding repeats to SONIC.\n");
    
  /* assuming vanilla RepeatMasker .out files here -- concatenated */
  while (fscanf(reps_file, "%s", sw_score) > 0){

    if (!strcmp(sw_score, "SW")){
      return_value_char = fgets(skip_line, MAX_LENGTH, reps_file); // skip the first header line
      return_value_char = fgets(skip_line, MAX_LENGTH, reps_file); // skip the second header line
      return_value_char = fgets(skip_line, MAX_LENGTH, reps_file); // skip the empty line
      continue;
    }
      
    /* data starts here */

    return_value = fscanf(reps_file, "%f%f%f%s%d%d%s%s%s%s%s%s%s%d\n", &perc_div, &perc_del, &perc_ins, chromosome,
			  &start, &end, chrom_left, strand, repeat_type, repeat_class, repeat_start_string, repeat_end_string, repeat_left, &repeat_id);

    /* write the necessary ones */
    chrom_name_length = strlen(chromosome);
    return_value = gzwrite(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzwrite(sonic_file, chromosome, chrom_name_length);
    return_value = gzwrite(sonic_file, &start, sizeof(start));
    return_value = gzwrite(sonic_file, &end, sizeof(end));
    return_value = gzwrite(sonic_file, strand, 1);
    repeat_type_length = strlen(repeat_type);
    repeat_class_length = strlen(repeat_class);
    return_value = gzwrite(sonic_file, &repeat_type_length, sizeof(repeat_type_length));
    return_value = gzwrite(sonic_file, repeat_type, repeat_type_length);
    return_value = gzwrite(sonic_file, &repeat_class_length, sizeof(repeat_class_length));
    return_value = gzwrite(sonic_file, repeat_class, repeat_class_length);
    if (strand[0] == '+'){
      repeat_start = atoi(repeat_start_string);
      repeat_end   = atoi(repeat_end_string);
    }
    else{
      repeat_start = atoi(repeat_left);
      repeat_end   = atoi(repeat_end_string);      
    }
    
    return_value = gzwrite(sonic_file, &repeat_start, sizeof(repeat_start));
    return_value = gzwrite(sonic_file, &repeat_end, sizeof(repeat_end));
    
  }
  
  fclose(reps_file);
  gzclose(sonic_file);

  fprintf( stderr, "SONIC file %s is ready.\n", sonic);
  return RETURN_SUCCESS;
  
}


  /* 
     TODO:  implement this. Replace write/read
  */

int load_sonic(char *sonic){

  gzFile sonic_file;
  int sonic_magic;
  int return_value;
  char *return_value_char;
  int line_count;
  int num_lines_read;
  char chromosome[MAX_LENGTH];
  int chrom_name_length;
  int start, end;

  /* 
     RepeatMasker out file entries
  */


  char strand[2];
  char repeat_type[MAX_LENGTH];
  char repeat_class[MAX_LENGTH];
  int repeat_type_length;
  int repeat_class_length;
  int repeat_start, repeat_end;
  
  fprintf (stderr, "Loading SONIC file..\n");
		     
  sonic_file = sonic_fopen_gz(sonic, "r");

  
  return_value = gzread(sonic_file, &sonic_magic, sizeof(sonic_magic));

  if (return_value == 0){
    fprintf(stderr, "Cannot read the SONIC file.\n");
    return EXIT_FILE_OPEN_ERROR;
  }

  if (sonic_magic != SONIC_MAGIC){
    fprintf(stderr, "Invalid SONIC file (%s). Please use a correctly created SONIC file.\n", sonic);
    return EXIT_SONIC;
  }


  return_value = gzread(sonic_file, &line_count, sizeof(line_count));
  fprintf( stderr, "Loading %d gap intervals.\n", line_count);
  num_lines_read = 0;
  
  while (num_lines_read != line_count){
    chrom_name_length = strlen(chromosome);
    return_value = gzread(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzread(sonic_file, chromosome, chrom_name_length);
    chromosome[chrom_name_length] = 0;
    return_value = gzread(sonic_file, &start, sizeof(start));
    return_value = gzread(sonic_file, &end, sizeof(end));
    /* fprintf(stderr, "[gaps]\t%s\t%d\t%d\n", chromosome, start, end); */
    num_lines_read++;
  }
  


  return_value = gzread(sonic_file, &line_count, sizeof(line_count));
  fprintf( stderr, "Loading %d duplication intervals.\n", line_count);
  num_lines_read = 0;

  while (num_lines_read != line_count){
    return_value = gzread(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzread(sonic_file, chromosome, chrom_name_length);
    chromosome[chrom_name_length] = 0;
    return_value = gzread(sonic_file, &start, sizeof(start));
    return_value = gzread(sonic_file, &end, sizeof(end));
    /* fprintf(stderr, "[dups]\t%s\t%d\t%d\n", chromosome, start, end); */
    num_lines_read++;
  }

  
  fprintf( stderr, "Loading repeats.\n");
  
  while (1){
    if (gzeof(sonic_file))
      break;

    chrom_name_length = strlen(chromosome);
    return_value = gzread(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
    return_value = gzread(sonic_file, chromosome, chrom_name_length);
    chromosome[chrom_name_length] = 0;
    return_value = gzread(sonic_file, &start, sizeof(start));
    return_value = gzread(sonic_file, &end, sizeof(end));
    return_value = gzread(sonic_file, strand, 1);
    strand[1] = 0;
    return_value = gzread(sonic_file, &repeat_type_length, sizeof(repeat_type_length));
    return_value = gzread(sonic_file, repeat_type, repeat_type_length);
    return_value = gzread(sonic_file, &repeat_class_length, sizeof(repeat_class_length));
    return_value = gzread(sonic_file, repeat_class, repeat_class_length);
    repeat_type[repeat_type_length] = 0;
    repeat_class[repeat_class_length] = 0;
    return_value = gzread(sonic_file, &repeat_start, sizeof(repeat_start));
    return_value = gzread(sonic_file, &repeat_end, sizeof(repeat_end));

    /* fprintf(stderr, "[reps]\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n", chromosome, start, end, strand, repeat_type, repeat_class, repeat_start, repeat_end); */

  }
  
  

  gzclose(sonic_file);
  return RETURN_SUCCESS;

  
}

FILE* sonic_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		fprintf( stderr, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		exit(EXIT_FILE_OPEN_ERROR);

	}
	return file;
}

gzFile sonic_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	gzFile file;
	char err[500];

	file = gzopen( path, mode);  
	if( !file)
	{
	        fprintf( stderr, "[TARDIS INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
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
		if (feof(bed_file))
			break;
		if (line[0] != 0)
			number_of_lines++;
	}

	return number_of_lines;
}
