#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M jggatter@mit.edu
#############################################

module add jellyfish/2.0.0

: '
if [ -f threemer_count_jellyfish.sh.e* ]; then
	rm threemer_count_jellyfish.sh.e*
	rm threemer_count_jellyfish.pe*
	rm threemer_count_jellyfish.o*
	rm threemer_count_jellyfish.sh.po*
fi
'

printf "Counting 3mers across genome...\n"
jellyfish count -m 3 -s 3300M -t 32 hg38_chromosomes.fa

printf "Making histogram of counts\n"
jellyfish histo mer_counts.jf

printf "Dumping counts...\n"
jellyfish dump -c mer_counts.jf > mer_counts_dump.fa

printf "DONE"