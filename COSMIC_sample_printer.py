import csv

SAMPLE_NAME=4
PRIMARY_HISTOLOGY=11
MUTATION_CDS=17
GRCH=22
MUTATION_GENOME_POSITION=23
MUTATION_STRAND=24
MUTATION_SOMATIC_STATUS=29
PUBMED_PMID=30

with open("V85_38_MUTANT.csv") as csvfile:
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
			#for x in range(0, 35): print(row[x], end=" ")
			print(row[PUBMED_PMID], end="")
			print("\n")

		#print("DEBUG:", row[SAMPLE_NAME], row[MUTATION_CDS], row[MUTATION_STRAND], row[MUTATION_GENOME_POSITION])
print("DONE.")
