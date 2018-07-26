#create_spectra.sh
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#mutpos_manual.py author: linakim [at] mit.edu
#July 26th, 2018

if [ $# -lt 3 ] || [ $# -gt 3 ]; then
	printf "ERROR: three arguments must be specified\n"
	printf "Format as: ./create_spectra.sh [samples_contexted_dir] [ref_count.txt] [output figures_dir]\n"
	printf "Example: ./create_spectra.sh ./samples_contexted/ ../ref/hg38_chromosomes_ref_counts.txt ./samples_contexted_figures\n"
	exit 1
fi

samples_contexted_dir=$1 #./samples_contexted/
ref_count=$2 #../ref/hg38_chromosomes_ref_counts.txt
figures_dir=$3 #./samples_contexted_figures
if [ ! -d "$figures_dir" ]; then 
	mkdir $figures_dir 
fi

samples_contexted=($samples_contexted_dir/*.csv)
for sample in ${samples_contexted[@]}; do

	sample_name=${sample#${samples_contexted_dir}}
	sample_name=${sample_name%.csv}	
	
	printf "Counting lines in $sample...\n"
	line_count=$(wc -l $sample | cut -d ' ' -f1)
	if [ $line_count -lt 51 ]; then 
		printf "	Skipping ${sample_name} due to line count of $line_count\n"
		continue 
	fi
	
	entry_count=$(( $line_count - 1 ))
	printf "$entry_count entries detected. Creating figure for ${sample_name}...\n"
	python mutpos_manual.py -k $ref_count -i $sample -s ${figures_dir}/${sample_name}
	printf "	Finished ${sample_name}.\n"
	#-k needs to be something else, the 3mer/trinucleotide frequencies within hg38

done