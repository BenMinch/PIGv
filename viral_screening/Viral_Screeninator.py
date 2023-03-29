###Diamond Blast Batch Script####
import sys,argparse,os,subprocess,re
import pandas as pd
import numpy as np
##flags
parser = argparse.ArgumentParser(description='''Script to use Diamond Blast in batch mode''')
parser.add_argument('-i', '--input', help='''Input folder of raw or trimmed reads''', required=True)

args = parser.parse_args()
input_folder = args.input
dataframe=pd.DataFrame(columns=['file','queries_aligned','Chance_Of_Finding_Genome'])
#make an index column with a number of rows equal to the number of files in the folder
dataframe['index']=range(1,len(os.listdir(input_folder))+1)
i=0
for file in os.listdir(input_folder):
    i+=1
    path= os.path.join(input_folder, file)
    print('Running Diamond Blast on '+file)
    cmd1='seqtk seq -A '+path+' > '+file+'.fasta'
    subprocess.call(cmd1, shell=True)
    cmd2='./diamond blastx -d MCP -q '+file+'.fasta --log'
    subprocess.call(cmd2, shell=True)
    cmd3='rm '+file+'.fasta'
    subprocess.call(cmd3, shell=True)
    grep= 'grep "queries aligned" diamond.log > '+file+'.txt'
    subprocess.call(grep, shell=True)
    cmd4='rm diamond.log'
    subprocess.call(cmd4, shell=True)
    #create a dataframe

    dataframe['file'][i-1]=file
    with open(file+'.txt') as f:
        split=f.read().split()[0]
    #make split numeric
    split=int(split)
    dataframe['queries_aligned'][i-1]=split
    #remove the file
    cmd5='rm '+file+'.txt'
    subprocess.call(cmd5, shell=True)
    #write the dataframe to a csv
dataframe.to_csv('viral_screening.csv', index=False)
    #run the r script
rscript= 'Rscript scripts/screening.r viral_screening.csv'
subprocess.call(rscript, shell=True)
