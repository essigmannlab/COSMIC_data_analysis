from pyfaidx import Fasta

key_names = []
for i in range(1,23):
	key_names.append("chr"+str(i))
key_names.append("chrX")
key_names.append("chrY")

#chromosomes = Fasta('hg38.fa')
#print(chromosomes.keys())
#print("Writing chromosomes")
'''
output = open('hg38_chromosomes.fa', 'a')
for key in key_names:
	print("Writing " + key)
	for line in chromosomes[key]:
		output.write(str(line)+"\n")
	output.write("\n")
output.close()
'''
print("hello")
new_file = Fasta('hg38_chromosomes.fa')
print(new_file.keys())
print("DONE")