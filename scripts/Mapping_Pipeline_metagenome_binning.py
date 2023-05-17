import os, sys, re, shlex, subprocess
import pandas as pd
#Filters out the internal standards
#cores = sys.argv[1]
input_dir = sys.argv[1]
outdir = sys.argv[2]
assembled_contigs= sys.argv[3]
#names of fasta and SRR run
reflist= sys.argv[4]


#read in the list of reflist

data= pd.read_csv(reflist)
#cycle through the reflist
os.mkdir(outdir)
for i in range(len(data)):
	fasta= data['File'][i]
	SRR= data['Run'][i]
	print(fasta)

#find fasta assembled contig file
	for n in os.listdir(assembled_contigs):
		if n==fasta:
			assembly= os.path.join(assembled_contigs, n)
		#cycle through the input directory
#check if SRR is in the input directory		
	for n in os.listdir(input_dir):
		if SRR in n:
			forward = os.path.join(input_dir, n)
			reverse= re.sub('_1', '_2',forward)
		#get out of loop if forward equals reverse		
			if forward==reverse:
				print('skip')
			else:
				outpath= os.path.join(outdir, fasta)
				outpath2= re.sub('.fa','.coverm', outpath)
				cmd = 'coverm contig --coupled '+ forward+ ' '+ reverse +' --reference '+ assembly+ ' -t 12 --methods metabat > '+ outpath2
				subprocess.call(cmd, shell=True)




