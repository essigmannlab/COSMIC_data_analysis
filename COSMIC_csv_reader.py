import csv
import re
import os

Names_positions = {} # Sample_name : Positions_substitutions{}
# { Sample Name : {Position : {Substitution : Count}} }
def write_dictionary(sample_name, chromosome, position, substitution, Names_positions=Names_positions):
	
	print("	DEBUG: Writing " + sample_name + " " + chromosome + " " + position + " " + substitution)
	position = str(chromosome) + ":" + str(position)
	Substitutions_counts = {} # strSubstitution_count : int Count
	Substitutions_counts[substitution] = 1
	Positions_substitutions = {} # strPosition : Substitutions_count{}
	Positions_substitutions[position] = Substitutions_counts
	
	if Names_positions.get(sample_name) != None:
		if Names_positions.get(sample_name).get(position) != None:
			if Names_positions.get(sample_name).get(position).get(substitution) != None: 
				Names_positions.get(sample_name).get(position)[substitution] += 1
			else: (Names_positions.get(sample_name).get(position))[substitution] = 1
		else: (Names_positions.get(sample_name))[position] = Substitutions_counts
	else: Names_positions[sample_name] = Positions_substitutions
	
	if Names_positions.get(sample_name).get(position).get(substitution) == None: print("	ERROR: " + sample_name + " " + position + " " + substitution + " was not printed!")
	return Names_positions;

SAMPLE_NAME=4
PRIMARY_HISTOLOGY=11
MUTATION_CDS=17
GRCH=22
MUTATION_GENOME_POSITION=23
MUTATION_STRAND=24
MUTATION_SOMATIC_STATUS=29

with open("V85_38_MUTANT.csv") as csvfile:
#with open("mutation_csv_reader_test.csv") as csvfile:
	
	reader = csv.reader(csvfile)
	for row in reader:
		
		#print("DEBUG: " + row[SAMPLE_NAME] + row[PRIMARY_HISTOLOGY] + row[MUTATION_CDS] + row[GRCH] + row[MUTATION_SOMATIC_STATUS])
		if row[0] == "GENE_NAME": continue
		if row[PRIMARY_HISTOLOGY] != "malignant_melanoma": continue
		if "del" in row[MUTATION_CDS] or "ins" in row[MUTATION_CDS]: continue
		if row[GRCH] == "null" or int(row[GRCH]) != 38: continue
		if "null" in row[MUTATION_STRAND]: continue
		if row[MUTATION_SOMATIC_STATUS] != "Confirmed somatic variant": continue

		Tandem_bases = ['AA>', 'AC>', 'AG>', 'AT>', 'CC>', 'CA>', 'CG>', 'CT>', 'GG>', 'GA>', 'GC>', 'GT>', 'TT>', 'TA>', 'TC>', 'TG>']
		Substitutions = ['A>G', 'A>C', 'A>T', 'C>G', 'C>A', 'C>T', 'G>C', 'G>A', 'G>T', 'T>G', 'T>A', 'T>C']
		position=row[MUTATION_GENOME_POSITION]
		strand=row[MUTATION_STRAND]
		contains_tandem = False
		for tandem in Tandem_bases:
			if tandem in row[MUTATION_CDS]:
				
				contains_tandem = True
				positions = re.split(':|-', position)
				cds_positions = re.split('G|C|T|A', row[MUTATION_CDS]) #maybe position = re.split... #maybe [_GCTA]
				substitutions = row[MUTATION_CDS].replace(cds_positions[0], "")
				if len(substitutions) != 5:
					print("	ERROR: Skipping " + row[SAMPLE_NAME] + " " + positions[0] +" "+ positions[1] + " and " + positions[2] + " " + substitutions)
					print("		DEBUG: " + substitutions + " in " + row[MUTATION_CDS] + " does not have a length of 5 characters.")
					break
				if substitutions[0] != substitutions[3]: 
					write_sub1 = strand+":"+substitutions[0]+">"+substitutions[3]
					Names_positions = write_dictionary(sample_name=row[SAMPLE_NAME], chromosome=positions[0], position=positions[1], substitution=write_sub1)
				if substitutions[1] != substitutions[4]:
					write_sub2 = strand+":"+substitutions[1]+">"+substitutions[4]
					Names_positions = write_dictionary(sample_name=row[SAMPLE_NAME], chromosome=positions[0], position=positions[2], substitution=write_sub2)
				break

		if contains_tandem == False:
			for sub in Substitutions:
				if sub in row[MUTATION_CDS]:

					positions = re.split(':|-', position)
					write_sub1 = strand+":"+sub
					Names_positions = write_dictionary(sample_name=row[SAMPLE_NAME], chromosome=positions[0], position=positions[1], substitution=write_sub1)
					break

dir_name = "./samples/"
if not os.path.exists(dir_name): os.makedirs(dir_name)
for sample_name in Names_positions:

	if "/" in sample_name: filename = sample_name.replace("/", "") 
	else: filename = sample_name

	path = os.path.join(dir_name, filename+".csv")
	output = open(path, 'w+')
	print("Writing to " + path)
	output.write("Position,Sample Name,Strand:Substitution,Count\n")
	
	for position in Names_positions.get(sample_name):
		#print(" DEBUG: " + position)
		for substitution in Names_positions.get(sample_name).get(position):
			#print(" DEBUG: " + substitution)
			count = Names_positions.get(sample_name).get(position).get(substitution)
			#print(" DEBUG: Writing " + str(position) + "," + str(sample_name) + "," + str(substitution) + "," + str(count))
			output.write(str(position) + "," + str(sample_name) + "," + str(substitution) + "," + str(count) + "\n")

print("Number of different samples: " + str(len(Names_positions)))
print("Done")