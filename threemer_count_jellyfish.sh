#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M jggatter@mit.edu
#############################################

module add jellyfish/2.0.0

#threemer_count_jellyfish.sh
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 26th, 2018

if [ $# -lt 1 ] || [ $# -gt 1 ]; then
	printf "ERROR: one argument must be specified\n"
	printf "Format as: ./create_spectra.sh [fasta file]\n"
	printf "Example: ./create_spectra.sh hg38_chromosomes.fa\n"
	exit 1
fi

fasta_file=$1

printf "Counting 3mers across genome...\n"
jellyfish count -m 3 -s 3300M -t 32 $fasta_file

printf "Making histogram of counts\n"
jellyfish histo mer_counts.jf

printf "Dumping counts...\n"
jellyfish dump -c mer_counts.jf > mer_counts_dump.fa

printf "DONE"