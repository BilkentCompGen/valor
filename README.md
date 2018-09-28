<<<<<<< HEAD
VALOR
======

VAriation with LOng Range information

Requirements
============

 * zlib   (http://www.zlib.net)
 * htslib (included; http://htslib.org/)
 * sonic  (included; https://github.com/calkan/sonic)

Fetching VALOR
===============

	git clone https://github.com/BilkentCompGen/valor.git --recursive

Compilation
===========

Type:
	
	make libs
	
	make
	
	cp valor /path/to/your/favorite/binaries (or sudo make install)


Running VALOR
==============

Required Parameters:

	-i, --input [BAM files]        : Input files in sorted BAM format.
	
	-o, --out   [output folder]    : Folder to put results and temporary files in
        
	-s, --sonic  [sonic file]      : Sonic file. Check: https://github.com/calkan/sonic.
        
	-f, --svs_to_find   [sv type]  : Comma separated list of SVs (i.e. INV,DUP,IDUP)

Optional Parameters:
      
	-l, --log_file [logfile name]: default is [output folder]/valor.log
	
Help Parameters:
        
	-v, --version                  : Print version and exit.
        
	-h, --help                     : Print help screen and exit.

	
Example:
	
	valor -i input.bam -s hg19.sonic -o folder_name -f DUP

SONIC file (annotations container)
==================================

SONIC files are available under https://github.com/BilkentCompGen/sonic-prebuilt/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)

 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.

 * GRCh38.sonic: SONIC file for the human reference genome build 38.
	
Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building a new SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

OUTPUT FORMAT
=============

in [OUTPUT_DIR]/predicted_inversions.bedpe

```bed
Chromosome-name BP1-start BP1-end Chromosome-name BP2-start BP2-end SV_TYPE 10XG-Support Read_Pair_Support Molecule_Depth
```
* Read pair Support: Number of read-pairs that support these breakpoints
* 10XG Support: Number of 10XG Molecule-pairs that support these breakpoints

Docker Usage
============

To build a valor Docker image

```
cd docker
docker build . -t valor:latest
```

Your image named "valor" should be ready. You can run valor using this image by

```
docker run --user=$UID -v /path/to/inputs:/input -v /path/to/outputdir:/output valor [args]
```

- ```[args]``` are usual arguments you would pass to valor executable. Be careful about mapping. You need to specify folders respective to container directory structure.
- You need to map host machine input and output directory to responding volume directories inside the container. These options are specified by '-v' argment.
- Docker works with root user by default. "--user" option saves your outputs.
=======
# sonic
Some Organism's Nucleotide Information Container

# Compilation

Library only:

	make
       
Standalone SONIC builder:

	make exe

# Standalone SONIC builder example:

	sonic --ref hg19.fasta --dups hg19.dups.bed --reps hg19_repeats.out --gaps hg19.gap.bed --make-sonic ucsc_hg19.sonic --info "UCSC_hg19"

	Duplications and gaps are expected in BED format. Repeats are in RepeatMasker .out format.

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
		GC% over the region, returned in [0-100] interval.
	
>>>>>>> dab878aca4f5fcae0149d3cd1a62e38eb3187cb3
