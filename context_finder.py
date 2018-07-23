import csv
import os
import sys
import signal
import glob
import re
from pyfaidx import Fasta

POSITION=0
SAMPLE_NAME=1
SUBSTITUTION=2
COUNT=3

success_counts={}
total_counts={}
for i in range(1,23): 
	success_counts[str(i)] = 0
	total_counts[str(i)] = 0
success_counts["X"] = 0
total_counts["X"] = 0
success_counts["Y"] = 0
total_counts["Y"] = 0
success_counts["M"] = 0
total_counts["M"] = 0
TOTAL="TOTAL"
total_counts[TOTAL] = 0
success_counts[TOTAL] = 0

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

fasta_filename = "hg38.fa"
print("Reading " +fasta_filename+ ", this may take some time...")
chromosomes = Fasta(fasta_filename)

print("Finished reading, cycling CSV files")
#samples = glob.glob("./samples_test_context/*.csv")
samples = glob.glob("./samples/*.csv")
dir_name = "./samples_contexted/"
if not os.path.exists(dir_name): os.makedirs(dir_name)
for sample_file in samples:
	
	print("In " + sample_file)	
	#sample_filename = sample_file.replace("./samples_test_context/","").replace(".csv","")
	sample_filename = sample_file.replace("./samples/","").replace(".csv","")
	contexted_sample_filename = sample_filename+"_contexted.csv"
	output_path = os.path.join(dir_name, contexted_sample_filename)
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
			
			start=position-2
			end=position+1
			strand_substitution = row[SUBSTITUTION].split(':')
			strand = strand_substitution[0]
			#print("		DEBUG: strand " + strand)
			base1 = (strand_substitution[1])[0]
			#print("		DEBUG: base1 " + base1 + " in " + strand_substitution[1])
			
			absolute_context = chromosomes[chromosome_entryname][start:end].seq
			#print("		DEBUG: absolute_context " + absolute_context)
			if strand == '+': 
				correct_base = base1
				true_context = absolute_context
			elif strand == '-':
				correct_base = reverse_complement(base1)
				#print("		DEBUG: reversed to correct base " + correct_base)
				true_context = reverse_complement(absolute_context)
				#print("		DEBUG: reversed to true context " + true_context)
			mutpos_abs = absolute_context[1]
			before_mutpos_abs = absolute_context[0]
			after_mutpos_abs = absolute_context[2]


			if  mutpos_abs.upper() != correct_base: 
				print("		ERROR: incorrect middle base "+before_mutpos_abs+mutpos_abs+after_mutpos_abs+" on the + strand when it should be "+correct_base)
				total_counts[chromosome]+=1
				total_counts[TOTAL]+=1
				#print("		DEBUG: total_counts: " + str(total_counts[TOTAL]))
				true_context = "INC_"+ true_context
			else: 
				success_counts[chromosome]+=1
				success_counts[TOTAL]+=1
				total_counts[chromosome]+=1
				total_counts[TOTAL]+=1
				print("		SUCCS: " + absolute_context +" contains "+ correct_base + " writing true context " + true_context)
				#print("		DEBUG: " + chromosome +" success_counts: " + str(success_counts[chromosome]) + " out of " + str(total_counts[chromosome]))
				#print("		DEBUG: Overall success_counts: " + str(success_counts[TOTAL]) + " out of " + str(total_counts[TOTAL]))
				
				output = open(output_path, 'a')
				output.write(str(row[POSITION]) + "," + true_context.upper() + "," + strand_substitution[1] + "," + str(row[COUNT]) + "\n")
				output.close()

	for chromosome in success_counts:
		if chromosome == TOTAL: continue
		if total_counts[chromosome] == 0: continue
		print(chromosome + " success rate: " + str(100*(success_counts[chromosome]/total_counts[chromosome])) + " percent from " + str(total_counts[chromosome]) + " entries.")
	print("Overall success rate: " + str(100*(success_counts[TOTAL]/total_counts[TOTAL])) + " percent from " + str(total_counts[TOTAL]) + " entries.")
	#print(chromosomes.keys())
	
print("DONE. Overall stats are shown above.")
#signal.pause()