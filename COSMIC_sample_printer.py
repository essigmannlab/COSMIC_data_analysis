import csv
import argparse

#COSMIC_sample_printer.py
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 26th, 2018

SAMPLE_NAME=4
PRIMARY_HISTOLOGY=11
MUTATION_CDS=17
GRCH=22
MUTATION_GENOME_POSITION=23
MUTATION_STRAND=24
MUTATION_SOMATIC_STATUS=29
PUBMED_PMID=30

parser = argparse.ArgumentParser(description="""A small portion of COSMIC_csv_reader.py that filters by certain 
											  parameters and prints the samples to console. Used primarily for 
											  searching for details of particular samples. You will need to hardcode
											  parameters by which you wish to search.""",
								 epilog="James was here")
parser.add_argument("-i", "--input",
					default="V85_38_MUTANT.csv"
					help="The input COSMIC .csv file")
args = parser.parse_args()

with open(args.input) as csvfile:
#with open("mutation_csv_reader_test.csv") as csvfile:
	
	reader = csv.reader(csvfile)
	for row in reader:
		
		if row[0] == "GENE_NAME": continue
		if row[PRIMARY_HISTOLOGY] != "malignant_melanoma": continue
		if "del" in row[MUTATION_CDS] or "ins" in row[MUTATION_CDS]: continue
		if row[GRCH] == "null" or int(row[GRCH]) != 38: continue
		if "null" in row[MUTATION_STRAND]: continue
		if row[MUTATION_SOMATIC_STATUS] != "Confirmed somatic variant": continue
		if str(row[SAMPLE_NAME]) == "2492704" or str(row[SAMPLE_NAME]) == "2492705" or str(row[SAMPLE_NAME]) == "2492706":
			print("FOUND " + row[SAMPLE_NAME])
			for x in range(0, 35): print(row[x], end=" ") #EVERYTHING
			#print(row[PUBMED_PMID], end="") #Prints only PUBMED_ID
			print("\n")
print("DONE.")
