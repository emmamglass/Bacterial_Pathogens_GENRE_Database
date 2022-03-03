outfile = open('job.txt','w')

with open('genomeids.txt', 'r') as f:
	genomeids = f.readlines()

with open('gram.txt', 'r') as f:
	grams = f.readlines()

new_genomeids = []
new_grams =[]
for ids in genomeids:
	new_genomeids.append(ids.strip())

#for gram in grams:
#	new_grams.append(gram.strip())


i = 0
for ids in new_genomeids:
	outp = 'srun python reconstruct_emg.py --input ' + str(new_genomeids[i]) + '.fa --type 1 --gram ' + str(grams[i])
	outfile.write(outp)
	i += 1
outfile.close()

