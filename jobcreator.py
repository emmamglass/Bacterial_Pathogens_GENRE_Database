import glob, os 
import natsort

outfile = open('job.txt','w')

genomeids = []
filelist = glob.glob("*.fa")
filelist = natsort.natsorted(filelist)#filelist.sort(key = lambda f: int(re.sub('\D', '', f)))

for file in filelist:
	num = os.path.splitext(file)[0]
	genomeids.append(num)

with open('gram.txt', 'r') as f:
	grams = f.readlines()

#new_genomeids = []
new_grams =[]
for gram in grams:
	new_grams.append(gram.strip())

i = 0
for ids in genomeids:
	outp = 'srun --ntasks 1 python reconstruct_emg.py --input ' + str(genomeids[i]) + '.fa --type 1 --gram ' + str(new_grams[i]) + ' &\n'
	outfile.write(outp)
	i += 1
outfile.close()
