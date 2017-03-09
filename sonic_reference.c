
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
    retval = fscanf( ref_index, "%s%d%*s%*s%*d", chrom_name, &chrom_len);
    if ( retval <= 0)
      return RETURN_ERROR;
    (*names)[i] = NULL;
    sonic_set_str( &((*names)[i]), chrom_name);
    (*lengths)[i] = chrom_len;    
  }
  
  return RETURN_SUCCESS;
}
