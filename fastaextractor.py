from Bio import SeqIO

for record in SeqIO.parse("output.fasta", "fasta"):
	species = record.description.split('|')[-1].strip("]").strip() #.split(".")[0]
	with open(f"{species}.fa", "a") as f:
		SeqIO.write(record, f, "fasta")