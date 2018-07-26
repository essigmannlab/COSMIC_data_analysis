import csv
import os
import sys
import signal
import glob
import re
import argparse
from pyfaidx import Fasta

#context_finder.py
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 26th, 2018

POSITION=0
SAMPLE_NAME=1
SUBSTITUTION=2
COUNT=3

success_counts={}
total_counts={}
for i in range(1,23): success_counts[str(i)] = total_counts[str(i)] = 0
success_counts["X"] = total_counts["X"] = 0
success_counts["Y"] = total_counts["Y"] = 0
success_counts["M"] = total_counts["M"] = 0
TOTAL="TOTAL"
total_counts[TOTAL] = success_counts[TOTAL] = 0

def reverse_complement(sequence):
	output=""
	if len(sequence) == 0: return output
	reverse_sequence = sequence[::-1].upper()
	for base in reverse_sequence:
		if base == 'G': output += 'C'
		elif base == 'C': output += 'G'
		elif base == 'A': output += 'T'
		elif base == 'T': output += 'A'
	return output
'''
def signal_handler(sig, frame, success_counts=success_counts):
	#print("\n")
	for chromosome in success_counts:
		if chromosome == TOTAL: continue
		if total_counts[chromosome] == 0: continue
		print(chromosome + " success rate: " + str(100*(success_counts[chromosome]/total_counts[chromosome])) + "percent from " + str(total_counts[chromosome]) + " entries.")
	print("Overall success rate: " + str(100*(success_counts[TOTAL]/total_counts[TOTAL])) + "percent from " + str(total_counts[TOTAL]) + " entries.")
	sys.exit(0) 
signal.signal(signal.SIGINT, signal_handler)
'''

parser = argparse.ArgumentParser(description="""Reads the .csv files outputted by COSMIC_csv_reader.py and uses 
												pyfaidx to locate the sequence context of the mutation position. 
												The program calculates and prints the percent of successful matches 
												of the original base to genomic position. As the program counts the 
												number of occurrences of a sequence context occurring at a certain 
												genomic position within an individual sample, a great amount of memory 
												is required and it is recommended that the script run with about 4 GB of 
												memory. The sequence contexts for each sample are outputted into .csv files 
												within a new directory. If a base does not match to the position, it is not 
												written to the outputted to its respective .csv file. The format of each 
												outputted file is different than those outputted by COSMIC_csv_reader.py to 
												allow create_spectra.sh (specifically mutpos_manual.py) to read and create 
												spectra.""")
parser.add_argument("-f", "--fasta",
					default="hg38_chromosomes.fa",
					help="The input .fasta file")
parser.add_argument("-i", "--input",
					default="./samples/",
					help="The path of the directory containing sample .csv files outputted by COSMIC_csv_reader.py.")
parser.add_argument("-o", "--output",
					default="./samples_contexted/",
					help="Directory for storing outputted .csv files that will be read by create_spectra.sh/create_uhc_dendrogram_heatmap.py")
args = parser.parse_args()

fasta_filename = args.fasta
print("Reading " + fasta_filename + ", this may take some time...")
chromosomes = Fasta(fasta_filename)
samples = glob.glob(args.input+"*.csv")
output_dirname = args.output
if not os.path.exists(output_dirname): os.makedirs(output_dirname)
print("Finished reading, cycling CSV files")
for sample_file in samples:
	
	print("In " + sample_file)	
	sample_filename = sample_file.replace("./samples/","").replace(".csv","")
	contexted_sample_filename = sample_filename+"_contexted.csv"

	output_path = os.path.join(output_dirname, contexted_sample_filename)
	output = open(output_path, 'w+')
	output.write("Position,Context,Substitution,Count\n")
	output.close()
	
	with open(sample_file) as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:

			if row[POSITION] == "Position": continue
			print("	DEBUG: " + row[POSITION])
			chromosome_position = row[POSITION].split(':')
			chromosome = chromosome_position[0]
			if chromosome == "23": chromosome = "X"
			if chromosome == "24": chromosome = "Y"
			if chromosome == "25": continue
			#if chromosome == "25": chromosome = "M"
			position = int(chromosome_position[1])
			chromosome_entryname = "chr"+chromosome
			#print(chromosomes[chromosome_entryname][position-2:position+1].fancy_name)
			
			strand_substitution = row[SUBSTITUTION].split(':')
			strand = strand_substitution[0]
			base1 = (strand_substitution[1])[0]
			start=position-2
			end=position+1
			absolute_context = chromosomes[chromosome_entryname][start:end].seq
			if strand == '+': 
				correct_base = base1
				true_context = absolute_context
			elif strand == '-':
				correct_base = reverse_complement(base1)
				true_context = reverse_complement(absolute_context)
			mutpos_abs = absolute_context[1]
			before_mutpos_abs = absolute_context[0]
			after_mutpos_abs = absolute_context[2]


			if  mutpos_abs.upper() != correct_base: 
				print("		ERROR: incorrect middle base "+before_mutpos_abs+mutpos_abs+after_mutpos_abs+" on the + strand when it should be "+correct_base)
				total_counts[chromosome]+=1
				total_counts[TOTAL]+=1
				true_context = "INC_"+ true_context
			else: 
				success_counts[chromosome]+=1
				success_counts[TOTAL]+=1
				total_counts[chromosome]+=1
				total_counts[TOTAL]+=1
				print("		SUCCS: " + absolute_context +" contains "+ correct_base + " writing true context " + true_context)
				output = open(output_path, 'a')
				output.write(str(row[POSITION]) + "," + true_context.upper() + "," + strand_substitution[1] + "," + str(row[COUNT]) + "\n")
				output.close()

	for chromosome in success_counts:
		if chromosome == TOTAL: continue
		if total_counts[chromosome] == 0 or total_counts[chromosome] == None: continue
		print(chromosome + " success rate: " + str(100*(success_counts[chromosome]/total_counts[chromosome])) + " percent from " + str(total_counts[chromosome]) + " entries.")
	print("Overall success rate: " + str(100*(success_counts[TOTAL]/total_counts[TOTAL])) + " percent from " + str(total_counts[TOTAL]) + " entries.")
	
print("DONE. Overall stats are shown above.")
#signal.pause()