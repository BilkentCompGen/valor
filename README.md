VALOR
======

Discover inversions using long range information

Requirements
============

 * zlib   (http://www.zlib.net)
 * htslib (included; http://htslib.org/)
 * sonic  (included; https://github.com/calkan/sonic)

Fetching VALOR
===============

	git clone git@github.com:BilkentCompGen/valor.git --recursive

Compilation
===========


	Change SV_TO_FIND variable in the make file to find different variants	

Type:
	
	make libs
	
	make
	
	cp valor /path/to/your/favorite/binaries (or sudo make install)


Running VALOR
==============

Required Parameters:

	-i, --input [BAM files]        : Input files in sorted BAM format.
	
	-o, --out   [output folder]    : Folder to put stuff in
        
	-s, --sonic  [sonic file]      : Sonic file. Check: https://github.com/calkan/sonic.
        
	-f, --svs_to_find   [sv type]: Among INV,DUP,IDUP. Multi SV discovery not implemented.

Optional Parameters:
      
	-l, --log_file [logfile name]: default is valor.log
	
Help Parameters:
        
	-v, --version                  : Print version and exit.
        
	-h, --help                     : Print help screen and exit.

	
Example:
	
	valor -i input.bam -s hg19.sonic -o folder_name -f DUP

SONIC file (annotations container)
==================================

SONIC files are available under https://github.com/BilkentCompGen/sonic-prebuilt/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * GRCh38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building the SONIC file
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
