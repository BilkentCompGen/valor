
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

void sonic_write_gc_profile(gzFile sonic_file, FILE *ref_file, int number_of_chromosomes, char **chromosome_names, int *chromosome_lengths)
{

  char ch;
  int char_count;
  int window_id;
  int gc;
  char gc_content;
  int chromosome_index;
  char this_chromosome[255];
  char line[MAX_LENGTH];
  int chrom_name_length;
  int return_value;
  int written_gc;
  int end_of_gc;

  end_of_gc = SONIC_END_OF_GC;
  chromosome_index = -1;
  window_id = 0;

  written_gc = 0;
  
  while (!feof(ref_file)){
    
    ch = fgetc(ref_file);
    if (ch == EOF)
      break;

    if (ch == '>'){
      fscanf(ref_file, "%s", this_chromosome);
      fgets(line, MAX_LENGTH, ref_file);


      /*
      if (written_gc != 0)
	fprintf (stderr, "len: %d, expected %d windows, written %d\n", chromosome_lengths[chromosome_index], (chromosome_lengths[chromosome_index] / (SONIC_GC_WINDOW)), written_gc);
      */
      
      chromosome_index = sonic_find_chromosome_index(chromosome_names, this_chromosome, number_of_chromosomes);
      if (chromosome_index != -1){
	/* fprintf (stderr, "GC Chromosome: [%d]\t%s (%s)\n", chromosome_index, this_chromosome, chromosome_names[chromosome_index]); */
	gzwrite(sonic_file, &chromosome_index, sizeof(chromosome_index));
      }
	/*
      chrom_name_length = strlen(this_chromosome);
      return_value = gzwrite(sonic_file, &chrom_name_length, sizeof(chrom_name_length));
      return_value = gzwrite(sonic_file, this_chromosome, chrom_name_length);      */
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
	/*
	fprintf(stderr, "charcnt: %d\twinid: %d\tgc: %d - %f", char_count, window_id, (int) gc_content ,  (100.0 * gc / SONIC_GC_WINDOW) );
	getc(stdin); */
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

    if (!strcmp(chromosome_names[i], this_chromosome))
      return i;
    
  }

  return -1;
}


void  sonic_read_gc_profile(gzFile sonic_file, sonic *this_sonic)
{
  int chromosome_index;
  char gc_content;
  int window_id;
}
