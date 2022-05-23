import glob, os 
import natsort

outfile = open('pfluxjob.txt','w')

filelist = glob.glob("*.sbml")
filelist = natsort.natsorted(filelist)#filelist.sort(key = lambda f: int(re.sub('\D', '', f)))


for file in filelist:
	outp = 'srun --ntasks 1 python quickflux.py --input ' + str(file) +  ' &\n'
	outfile.write(outp)
outfile.close()
