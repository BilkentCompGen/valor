
#include "sonic_reference.h"

int get_number_of_chromosomes(FILE *ref_index){
  return count_bed_lines(ref_index);
}

int get_chromosome_info(FILE *ref_index, int number_of_chromosomes, int **lengths, char ***names){

  int i;
  int retval;
  char chrom_name[255];
  int chrom_len;
  
  (*lengths) = (int *) sonic_get_mem(sizeof(int) * number_of_chromosomes);
  (*names) = (char **) sonic_get_mem(sizeof(char *) * number_of_chromosomes);

  for (i=0; i < number_of_chromosomes; i++){
    retval = fscanf( ref_index, "%s%d%*s%*s%*d\n", chrom_name, &chrom_len);
    if ( retval <= 0)
      return RETURN_ERROR;
    (*names)[i] = NULL;
    sonic_set_str( &((*names)[i]), chrom_name);
    (*lengths)[i] = chrom_len;
  }
  
  return RETURN_SUCCESS;
}

void sonic_write_gc_profile(gzFile sonic_file, FILE *ref_file, int number_of_chromosomes, char **chromosome_names)
{

  char ch;
  int char_count;
  int window_id;
  int gc;
  char gc_content;
  int chromosome_index;
  char this_chromosome[255];
  char line[MAX_LENGTH];
  int return_value;
  char *return_value_char;
  int written_gc;
  int end_of_gc;

  end_of_gc = SONIC_END_OF_GC;
  chromosome_index = -1;
  window_id = 0;

  written_gc = 0;
  gc = 0;
  char_count = 0;
  
  while (!feof(ref_file)){
    
    ch = fgetc(ref_file);
    if (ch == EOF)
      break;

    if (ch == '>'){
      return_value = fscanf(ref_file, "%s", this_chromosome);
      if (return_value == EOF)
	break;
      
      return_value_char = fgets(line, MAX_LENGTH, ref_file);
      if (return_value_char == NULL)
	break;
      
      
      chromosome_index = sonic_find_chromosome_index(chromosome_names, this_chromosome, number_of_chromosomes);
      if (chromosome_index != -1){
	gzwrite(sonic_file, &chromosome_index, sizeof(chromosome_index));
      }
      char_count = 1;
      gc = 0;
      written_gc = 0;
      window_id = 0;
    }

    else if (!isspace(ch)){
      
      ch = toupper(ch);
      char_count++;

      if (ch == 'G' || ch == 'C')
	gc++;

      if (char_count % SONIC_GC_WINDOW == 0){
	gc_content = (char) (100.0 * gc / SONIC_GC_WINDOW);
	written_gc++;
	window_id++;
	gc = 0;
	gzwrite(sonic_file, &gc_content, sizeof(gc_content));
      }
    }
    
  }

  gzwrite(sonic_file, &end_of_gc, sizeof(end_of_gc));
}

int sonic_find_chromosome_index(char **chromosome_names, char *this_chromosome, int number_of_chromosomes){
  int i;

  /* linear scan. Not really proud of it, but it shouldn't affect performance anyway */
  for (i = 0; i < number_of_chromosomes; i++){
    if (!strcmp(chromosome_names[i], this_chromosome)){
      return i;
    }
  }

  return -1;
}

int sonic_refind_chromosome_index(sonic *this_sonic, char *this_chromosome){
  int i;

  if (this_sonic->last_chromosome_index != -1 && !strcmp(this_sonic->chromosome_names[this_sonic->last_chromosome_index], this_chromosome))
    return this_sonic->last_chromosome_index;
  
  /* linear scan. Not really proud of it, but it shouldn't affect performance too much */
  for (i = 0; i < this_sonic->number_of_chromosomes; i++){
    if (!strcmp(this_sonic->chromosome_names[i], this_chromosome)){
      this_sonic->last_chromosome_index = i;
      return i;
    }
  }

  return -1;
}


void sonic_read_gc_profile(gzFile sonic_file, sonic *this_sonic)
{
  int chromosome_index;
  char gc_content;
  int window_id;
  int return_value;
  int number_of_gc_windows;
  
  
  while (!gzeof(sonic_file)){
    return_value = gzread(sonic_file, &chromosome_index, sizeof(chromosome_index));
    if (chromosome_index == SONIC_END_OF_GC || return_value == EOF)
      break;

    number_of_gc_windows = this_sonic->chromosome_lengths[chromosome_index] / (SONIC_GC_WINDOW);
    
    window_id = 0;
    while (window_id < number_of_gc_windows){
      return_value = gzread(sonic_file, &gc_content, sizeof(gc_content));
      this_sonic->chromosome_gc_profile[chromosome_index][window_id++] = gc_content;
    }
  }

}
