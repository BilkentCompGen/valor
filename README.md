VALOR
======

Discover inversions using long range information

Requirements
============

 * zlib   (http://www.zlib.net)
 * htslib (included; http://htslib.org/)
 * sonic  (included; https://github.com/calkan/sonic)

Fetching valor
===============

	git clone git@github.com:BilkentCompGen/valor.git --recursive

Compilation
===========

Type:

	make libs
	make
	cp valor /path/to/your/favorite/binaries


Running VALOR
==============

	valor input.bam snc.sonic

	More command-line to be implemented ( Those can be changed from valorconfig.h for now)



SONIC file (annotations container)
==================================

SONIC files are available under aux/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building the SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

However, you can still build the SONIC file using TARDIS:

	tardis --ref human_g1k_v37.fasta --make-sonic human_g1k_v37.sonic \
		--dups human_g1k_v37.segmental_duplications.bed \
		--gaps human_g1k_v37.assembly_gaps.bed \
		--reps human_g1k_v37.repeatmasker.out 

	

