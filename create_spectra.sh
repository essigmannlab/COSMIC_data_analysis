samples_contexted=(./samples_contexted/*.csv)
for sample in ${samples_contexted[@]}; do

	sample_name=${sample#./samples_contexted/}
	sample_name=${sample_name%.csv}	
	
	printf "Counting lines in $sample...\n"
	line_count=$(wc -l $sample | cut -d ' ' -f1)
	if [ $line_count -lt 51 ]; then 
		printf "	Skipping ${sample_name} due to line count of $line_count\n"
		continue 
	fi
	
	entry_count=$(( $line_count - 1 ))
	printf "$entry_count entries detected. Creating figure for ${sample_name}...\n"
	python mutpos_manual.py -k ../ref/hg38_chromosomes_ref_counts.txt -i $sample -s ./samples_contexted_figures/${sample_name}
	printf "	Finished ${sample_name}.\n"
	#-k needs to be something else, the 3mer/trinucleotide frequencies within hg38

done