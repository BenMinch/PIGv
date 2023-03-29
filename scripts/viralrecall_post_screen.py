#viral recall batch genome screening
import pandas as pd
import numpy as np
import os, sys, re, shlex, subprocess
#a csv file with list of genomes (can be an ncldv output file)
masterfile= sys.argv[1]
#input folder of all summary files
input_folder= sys.argv[2]
#read in the masterfile
data= pd.read_csv(masterfile, sep='\t')
data['Sum_score']=''
data['Total_contigs']=''
data['Percent_negative']=''
data['Virus']=''
data['Contamination']=''
for i in range(len(data)):
    file= re.sub('.fa','.summary.tsv', data['genome'][i])
    for n in os.listdir(input_folder):
        if file == n:
            data2= pd.read_csv(os.path.join(input_folder, n), sep='\t')
            #check sum of a column in data2
            if data2['score'].sum() > 1:
                data['Virus'][i]= 'yes'
            else:
                data['Virus'][i]= 'no'
            data['Sum_score'][i]= data2['score'].sum()
            #count number of negative values in data2
            data['Contamination'][i]= len(data2[data2['score'] < 0])
            data['Total_contigs'][i]= len(data2)
            data['Percent_negative'][i]=data['Contamination'][i]/data['Total_contigs'][i]
# Optional Viral screening based on Aylward paper
data['PolB']=''
data['Other_Markers']=0
for i in range(len(data)):
    if data['GVOGm0054'][i] > 0:
        data['PolB'][i]= 'yes'
    else:
        data['PolB'][i]= 'no'
    if data['GVOGm0003'][i] > 0:
        data['Other_Markers'][i]= data['Other_Markers'][i] + 1
    if data['GVOGm0760'][i] > 0:
        data['Other_Markers'][i]= data['Other_Markers'][i] + 1
    if data['GVOGm0013'][i] >0:
        data['Other_Markers'][i]= data['Other_Markers'][i] + 1
    if data['GVOGm0890'][i] > 0:
        data['Other_Markers'][i]= data['Other_Markers'][i] + 1


#Write the data to csv
data.to_csv('viralrecall_batch_summary.tsv', sep='\t', index=False)


