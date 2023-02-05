import pysam
all_fasta = pysam.FastaFile("all.human.protein.fasta")

all_fasta.references[0]
all_fasta.fetch(all_fasta.references[0])

for tmp_protein_name in all_fasta.references:
	tmp_protein_seq = all_fasta.fetch(tmp_protein_name)
	# just use protein ID as file name, good for shell process
	tmp_fasta_name = "individual/" + tmp_protein_name
	# safe open and close file
	with open(tmp_fasta_name, mode="w") as fout:
		fout.write('>' + tmp_protein_name + '\n')
		fout.write(tmp_protein_seq + '\n')
	# break

# fout = open(tmp_fasta_name, mode='w')
# fout.close()

