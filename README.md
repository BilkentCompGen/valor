VALOR
======

VAriation with LOng Range information

*This repository includes code for VALOR2. The original version 1 is deprecated but also available in the rasekh-etal directory for archival purposes.*

Citation
============
	VALOR2: characterization of large-scale structural variants using linked-reads. 
	Fatih Karaoglanoglu, Camir Ricketts, Ezgi Ebren, Marzieh Eslami Rasekh,
	Iman Hajirasouliha, Can Alkan. Genome Biology, 21 (1): 72, 2020.

Requirements
============

 * gcc > 5.4
 * zlib   (http://www.zlib.net)
 * htslib (included; http://htslib.org/)
 * sonic  (included; https://github.com/calkan/sonic)

How to install VALOR from source
====================

	git clone https://github.com/BilkentCompGen/valor.git --recursive
	cd valor
	make libs
	make
	cp valor /path/to/your/favorite/binaries
	
Running VALOR
==============

Required Parameters:

	-i, --input [BAM files]        : Input files in sorted BAM format.
	
	-o, --out   [output prefix]    : Prefix for the output files
        
	-s, --sonic  [sonic file]      : Sonic file. Check: https://github.com/calkan/sonic.
        
	-f, --svs_to_find   [sv type]  : Comma separated list of SVs (i.e. INV,DUP,IDUP,TRA,ITRA,DEL)

Optional Parameters:
      
	-l, --log_file [logfile name]: default is [output prefix]-valor.log
	-t, --threads [thread count]: sets number of threads to be used (default is 1)
    -p, --ploidy  [Number of chromosome sets]      : Default is 2
    -y, --single-copy-chr  [Chromosome name]      : sets copy number of a chromosome to 1 (example -y X -y Y for male humans).	
Help Parameters:
        
	-v, --version                  : Print version and exit.
        
	-h, --help                     : Print help screen and exit.

	
Example:
	
	valor -i input.bam -s hg19.sonic -o output_prefix -f DUP,IDUP,DEL

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

in [output prefix]-predicted_inversions.bedpe

```bed
Chromosome-name BP1-start BP1-end Chromosome-name BP2-start BP2-end SV_TYPE 10XG-Support Read_Pair_Support Molecule_Depth
```
* Read pair Support: Number of read-pairs that support these breakpoints
* 10XG Support: Number of 10XG Molecule-pairs that support these breakpoints

BEDPE format table
==================
| chr1       	| start1         	| end1       	| chr2       	| start2       	| end2         	| type          	| Support info 	|
|------------	|----------------	|------------	|------------	|--------------	|--------------	|---------------	|--------------	|
| chr        	| bp1-start      	| bp1-end    	| chr        	| bp2-start    	| bp-end       	| inversion     	| ...          	|
| chr        	| molecule-start 	| start      	| chr        	| end          	| molecule-end 	| deletion      	|              	|
| chr        	| source-start   	| source-end 	| chr        	| target-start 	| target-end   	| duplication   	|              	|
| source-chr 	| source-start   	| source-end 	| target-chr 	| target-start 	| target-end   	| translocation 	|              	|


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
