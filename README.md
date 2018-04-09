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


	Change SV_TO_FIND variable in the make file to find different variants #TODO change this to cmdline variable	

Type:
	make libs
	make
	cp valor /path/to/your/favorite/binaries (or sudo make install)


Running VALOR
==============

	valor input.bam snc.sonic

	More command-line to be implemented ( Those can be changed from valorconfig.h for now)



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
Chromosome-name BP1-start BP1-end Chromosome-name BP2-start BP2-end SV_TYPE 10XG-Support Read_Pair_Support
```
* Read pair Support: Number of read-pairs that support these breakpoints
* 10XG Support: Number of 10XG Molecule-pairs that support these breakpoints
