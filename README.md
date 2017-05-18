# sonic
Some Organism's Nucleotide Information Container


# Building a SONIC file

	int sonic_build(char *ref_genome, char *gaps, char *reps, char *dups, char *info, char *sonic);

	ref_genome: path to the reference genome file (FASTA) [input]. Also requires the ref_genome.fai file in the same directory (samtools faidx).
	gaps: path to the gaps file (BED) [input].
	reps: path to the repeat annotation file (UCSC/RepeatMasker .out) [input].
	dups: path to the segmental duplications annotation file (BED) [input].
	info: information string to annotate the SONIC file [input].
	sonic: path to the SONIC file [output].
	
	RETURN VALUES:
		1: success
		5: file open error
		7: annotation files error
		

# Loading a SONIC file

	sonic *sonic_load(char *sonic_file_name);

	sonic_file_name: path to the SONIC file [input].

	RETURN VALUE:
		a SONIC.

# General-purpose interval search

	sonic_interval *sonic_intersect(sonic *sonic, char *chromosome, int start, int end, sonic_interval_type interval_type);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).
	interval_type:
		SONIC_GAP: search gap annotation.
		SONIC_DUP: search segmental duplication annotation.
		SONIC_REP: search repeats annotation.

	RETURN VALUE:
		pointer to a SONIC interval data structure on success.
		NULL if not found.

# Check if interval hits a satellite region

	int sonic_is_satellite(sonic *sonic, char *chromosome, int start, int end);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).

	RETURN VALUE:
		1 if the interval is in a satellite region.
		0 if it is not.

# Check if interval hits a gap

	int sonic_is_gap(sonic *sonic, char *chromosome, int start, int end);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).

	RETURN VALUE:
		1 if the interval is in a gap region.
		0 if it is not.

# Check if interval hits a segmental duplication region

	int sonic_is_segmental_duplication(sonic *sonic, char *chromosome, int start, int end);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).

	RETURN VALUE:
		1 if the interval is in a segmental duplication region.
		0 if it is not.

# Check if interval hits a mobile element

	sonic_repeat *sonic_is_mobile_element(sonic *sonic, char *chromosome, int start, int end, char *mei_string);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).
	mei_string: a colon-seperated of MEI keywords, i.e. "Alu:L1:HERV".

	RETURN VALUE:
		pointer to a SONIC repeat data structure if the interval hits a mobile element.
		NULL if it is not.

# Get GC content over a region (as pre-calculated in non-overlapping 100 bp intervals).

	float sonic_get_gc_content(sonic *sonic, char *chromosome, int start, int end);

	sonic: Loaded SONIC.
	chromosome: name of the chromosome of the interval of interest.
	start: start coordinate of the interval of interest (inclusive -- BED-like).
	end: end coordinate of the interval of interest (exclusive -- BED-like).

	RETURN VALUE:
		GC% over the region.
	
