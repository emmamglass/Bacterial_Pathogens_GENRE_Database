from Bio import SeqIO

species_orig = 'none'

#for every record in the sequence record
for record in SeqIO.parse("rrna_features.fasta", "fasta"):

	#grab the species name of the record
	species = str(record.description.split('|')[-1].strip("]").strip())

	#if we are dealing with a 16S sequence record
	if '16S' in record.description:

		#AND if we don't already have a 16S sequence from that given species
		if species != species_orig:
			print(species)
			print(species_orig)
			#Add that record to the 16S fasta record
			with open("16s_only.fasta", "a") as f:
				SeqIO.write(record, f, "fasta")

	species_orig = species

count = 0
for record in SeqIO.parse("16s_only.fasta", 'fasta'):
	count +=1
	print(count)