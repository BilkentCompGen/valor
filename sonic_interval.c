
#include "sonic_interval.h"
#include "sonic.h"

int sonic_insert(char *chromosome, int start, int end, sonic_interval_type sonic_type, sonic_repeat *repeat_item){

}


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
