import os,argparse,re,subprocess
import pandas as pd
import numpy as np


#####NCLDViWrap

## Arguments
parser = argparse.ArgumentParser(description='''A wrapper for mining NCLDV genomes from metagenomic assemblies''')
parser.add_argument('-i', '--input', help='Input metagenomic assembly', required=True)
parser.add_argument('-o', '--output', help='Output directory (it will make one with this name)', required=True)
parser.add_argument('-t', '--threads', help='Number of threads to use', required=True)
parser.add_argument('-cov', '--coverage', help='Input coverage file from mapping reads to assemblies', required=True)
parser.add_argument('-annot', '--annotation', help='Add this flag if you want annotations for genomes (True)', required=False)

args = parser.parse_args()
input=args.input
threads=args.threads
coverage=args.coverage
output=args.output
annotation=args.annotation
## Make output directory
os.mkdir(output)

## Markersearch on the contigs
print('Running NCLDV Markersearch on Contigs')
os.mkdir(output+'/00_Contig_Markersearch')
prodigal= 'python scripts/prodigal_launcher_single.py '+input+' '+output+'/00_Contig_Markersearch'
subprocess.call(prodigal, shell=True)
markersearch='python scripts/ncldv_markersearch_lowmcp.py -i '+output+'/00_Contig_Markersearch -a -m GVOGm0003 -t '+threads+' -n mcp'
subprocess.call(markersearch, shell=True)
move_mcp= 'mv mcp* '+output+'/00_Contig_Markersearch'
subprocess.call(move_mcp,shell=True)
## Bin the contigs
print('Binning Contigs')
os.mkdir(output+'/01_Metabat_bins')
command= 'metabat2 -i '+input+' -a '+ coverage+ ' -o '+output+'/01_Metabat_bins/bin -t '+threads+ ' -m 5000 -s 50000'
subprocess.call(command, shell=True)
print('Done')
## Run Prodigal on the bins
print('Running Prodigal')
prodigal_input= output+'/01_Metabat_bins'
os.mkdir(output+'/02_Prodigal_bins')
command= 'python scripts/prodigal_launcher.py '+prodigal_input+' '+output+'/02_Prodigal_bins'
subprocess.call(command, shell=True)
print('Done')
###### Run Markersearch on the Prodigals of the bins (remember to change markersearch)#####
print('Running NCLDV Markersearch')
os.mkdir(output+'/03_NCLDV_Markersearch')
os.mkdir(output+'/03_NCLDV_Markersearch/Predicted_proteins')
move= 'cp '+output+'/02_Prodigal_bins/*.faa '+output+'/03_NCLDV_Markersearch/Predicted_proteins'
subprocess.call(move, shell=True)

search= 'python scripts/ncldv_markersearch.py -i '+output+'/03_NCLDV_Markersearch/Predicted_proteins -a -m  GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013,GVOGm0022,GVOGm0003 -t '+threads+' -n '+output+'/03_NCLDV_Markersearch/markersearch'
subprocess.call(search, shell=True)
print('Done')
#####Sort through the tsv files and make a list of the bins with at least 1 hit#####
ncldv_out= pd.read_csv(output+'/03_NCLDV_Markersearch/markersearch.table.tsv', sep='\t')

#create a column called total that counts the number of nonzero values in each row
ncldv_out['total']=ncldv_out.iloc[:,1:].astype(bool).sum(axis=1)
#filter to only include rows with at least 2 hits
ncldv_out=ncldv_out[ncldv_out['total']>=2]
#write the genome names to a txt file
ncldv_out2=ncldv_out.iloc[:,0]
ncldv_out.to_csv(output+'/03_NCLDV_Markersearch/ncldv_summary.tsv', sep='\t', index=False)
ncldv_out2.to_csv(output+'/03_NCLDV_Markersearch/ncldv_bins.txt', sep='\t', index=False)
#extract the bins with at least 1 hit
os.mkdir(output+'/03_NCLDV_Markersearch/ncldv_good_bins')
extract= 'cat '+output+'/03_NCLDV_Markersearch/ncldv_bins.txt | xargs -I {} cp '+output+'/01_Metabat_bins/{} '+output+'/03_NCLDV_Markersearch/ncldv_good_bins'
subprocess.call(extract, shell=True)
## Run Viral recall on the bins with at least 1 hit
print('Running ViralRecall')
os.mkdir(output+'/04_Viral_recall')
recall= 'python scripts/viralrecall.py -i '+output+'/03_NCLDV_Markersearch/ncldv_good_bins -p '+output+'/04_Viral_recall/viral_recall_out -t '+threads+ ' -c -b'
subprocess.call(recall, shell=True)
print('Done')
## Run viral recall batch screen to screen for NCLDV
#get all the summary files
print('Screening bins')
os.mkdir(output+'/05_Viral_screening')
os.mkdir(output+'/05_Viral_screening/summary_files')
extract= 'python scripts/extract_anything.py '+output+'/04_Viral_recall/viral_recall_out '+output+'/05_Viral_screening/summary_files .summary.tsv'
subprocess.call(extract, shell=True)
post_screen= 'python scripts/viralrecall_post_screen.py '+ output+'/03_NCLDV_Markersearch/ncldv_summary.tsv'+ ' '+output+'/05_Viral_screening/summary_files'
subprocess.call(post_screen, shell=True)
#move the output file into the output directory
move= 'mv viralrecall_batch_summary.tsv '+output+'/05_Viral_screening'
subprocess.call(move, shell=True)
##move all the good genomes into a directory
os.mkdir(output+'/05_Viral_screening/ncldv_good_genomes')
viral_screen=pd.read_csv(output+'/05_Viral_screening/viralrecall_batch_summary.tsv', sep='\t')
viral_screen=viral_screen[viral_screen['Virus']=='yes']
viral_screen=viral_screen[viral_screen['Other_Markers']>=3]

viral_genomes=viral_screen.iloc[:,0]
viral_genomes.to_csv(output+'/05_Viral_screening/ncldv_good_genomes.txt', sep='\t', index=False)
extract= 'cat '+output+'/05_Viral_screening/ncldv_good_genomes.txt | xargs -I {} cp '+output+'/03_NCLDV_Markersearch/ncldv_good_bins/{} '+output+'/05_Viral_screening/ncldv_good_genomes'
subprocess.call(extract, shell=True)
print('Done')
####Remove contamination preprocessing####
os.mkdir(output+'/05_Viral_screening/clean_genomes')
#get list of files in ncldv_good_genomes
print('Decontaminating genomes')
files=os.listdir(output+'/05_Viral_screening/ncldv_good_genomes')
#make a variable with changed filenames to be .summary instead of .fa
new_files=[]
for file in files:
    new_files.append(file.replace('.fa','.summary.tsv'))
new_files=pd.DataFrame(new_files)
new_files.to_csv(output+'/05_Viral_screening/new_files.txt', sep='\t', index=False)
os.mkdir(output+'/05_Viral_screening/good_summary_files')
extract= 'cat '+output+'/05_Viral_screening/new_files.txt | xargs -I {} cp '+output+'/05_Viral_screening/summary_files/{} '+output+'/05_Viral_screening/good_summary_files'
subprocess.call(extract, shell=True)
#Parsing out the bad contigs in each genome
for file in os.listdir(output+'/05_Viral_screening/good_summary_files'):
    #read in the summary file if it ends with tsv
    if file.endswith('.tsv'):
        summary=pd.read_csv(output+'/05_Viral_screening/good_summary_files/'+file, sep='\t')
        summary=summary[summary['score']<=0]
        summary_removes= summary['replicon']
        summary_removes.to_csv(output+'/05_Viral_screening/good_summary_files/'+file+'.removes.txt', sep='\t', index=False)
####Removing the bad contigs from the genomes####
for file in os.listdir(output+'/05_Viral_screening/ncldv_good_genomes'):
    #read in the summary file if it ends with tsv
    if file.endswith('.fa'):
        #get the name of the file without the .fa
        name=file.replace('.fa','')
        #make a variable for the removes file
        removes=output+'/05_Viral_screening/good_summary_files/'+name+'.summary.tsv.removes.txt'
        #run the remove script
        remove= 'python scripts/Decontaminate_fastas.py '+output+'/05_Viral_screening/ncldv_good_genomes/'+file+' '+removes+' > '+output+'/05_Viral_screening/clean_genomes/'+file
        subprocess.call(remove, shell=True)
print('Done')
###### Gathering Final Results
os.mkdir(output+'/06_Final_Results')
move_dRep= 'mv '+output+'/05_Viral_screening/clean_genomes '+output+'/06_Final_Results'
subprocess.call(move_dRep, shell=True)
viral_screen.to_csv(output+'/06_Final_Results/viral_genome_statistics.tsv', sep='\t', index=False)
#### Clean up

remove1= 'rm '+ output+'/03_NCLDV_Markersearch/ncldv_bins.txt'
remove2= 'rm '+ output+'/03_NCLDV_Markersearch/ncldv_summary.tsv'
subprocess.call(remove1, shell=True)
subprocess.call(remove2, shell=True)
remove3= 'rm '+ output+'/05_Viral_screening/ncldv_good_genomes.txt'
remove4= 'rm '+ output+'/05_Viral_screening/new_files.txt'
subprocess.call(remove3, shell=True)
subprocess.call(remove4, shell=True)
remove5= 'rm -r '+ output+'/05_Viral_screening/good_summary_files'
subprocess.call(remove5, shell=True)


### GVclass taxonomy
print('Running GVclass')
## rename all genomes from .fa to .fna
for file in os.listdir(output+'/06_Final_Results/clean_genomes'):
    #read in the summary file if it ends with tsv
    if file.endswith('.fa'):
        #get the name of the file without the .fa
        name=file.replace('.fa','')
        name=name.replace('.','_')
        #make a variable for the removes file
        rename= 'mv '+output+'/06_Final_Results/clean_genomes/'+file+' '+output+'/06_Final_Results/clean_genomes/'+name+'.fna'
        subprocess.call(rename, shell=True)

###Gv class
gvlass= 'snakemake -j 24 --use-conda --config querydir="'+output+'/06_Final_Results/clean_genomes"'
subprocess.call(gvlass, shell=True)
move_gvclass= 'mv '+output+'/06_Final_Results/clean_genomes/results/gvclass_out.tab '+output+'/06_Final_Results'
subprocess.call(move_gvclass, shell=True)

##more cleanup
remove6= 'rm *.txt'
subprocess.call(remove6, shell=True)
remove7= 'rm *.err'
subprocess.call(remove7, shell=True)
remove8= 'rm *.out'
subprocess.call(remove8, shell=True)

###Optional annotations

if annotation == 'True':
    print('Annotating genomes')
    os.mkdir(output+'/06_Final_Results/annotations')
    prodigal_input= output+'/06_Final_Results/clean_genomes'
    command= 'python scripts/prodigal_launcher.py '+prodigal_input+' '+output+'/06_Final_Results/annotations'
    subprocess.call(command, shell=True)
    for file in os.listdir(output+'/06_Final_Results/annotations'):
        if file.endswith('.faa'):
            path= output+'/06_Final_Results/annotations/'+file
            path2= output+'/06_Final_Results/annotations/'
            name= file.replace('.faa','')
            name2=os.path.join(path2, name)
            hmmscan= 'hmmscan -E 0.001 --tblout '+name2+ '.hmmtable hmm/gvog.hmm '+ path
            subprocess.call(hmmscan, shell=True)
    ##parsing the outputs
    for file in os.listdir(output+'/06_Final_Results/annotations'):
        if file.endswith('.hmmtable'):
            path= output+'/06_Final_Results/annotations/'+file
            path2= output+'/06_Final_Results/annotations/'
            name= file.replace('.hmmtable','')
            name2=os.path.join(path2, name)
            parse= "awk '!x[$3]++' "+path+' > '+name2+'.parsed.tsv'
            subprocess.call(parse, shell=True)
    for file in os.listdir(output+'/06_Final_Results/annotations'):
        if file.endswith('.parsed.tsv'):
            path= output+'/06_Final_Results/annotations/'+file
            path2= output+'/06_Final_Results/annotations/'
            name= file.replace('.parsed.tsv','')
            name2=os.path.join(path2, name)
            parsinator= 'python scripts/hmm_parsinator.py '+path+' '+name2+'.parsinator.csv'
            subprocess.call(parsinator, shell=True)
    ##annotate the hmm
    for file in os.listdir(output+'/06_Final_Results/annotations'):
        if file.endswith('.parsinator.csv'):
            path= output+'/06_Final_Results/annotations/'+file
            path2= output+'/06_Final_Results/annotations/'
            name= file.replace('.parsinator.csv','')
            name2=os.path.join(path2, name)
            annotate= 'python scripts/hmm_annotationator.py -i '+path+' -o '+name2+'.annotated.csv -r hmm/gvog_annotation.tsv'
            subprocess.call(annotate, shell=True)
remove= 'rm '+output+'/06_Final_Results/annotations/*.hmmtable'
subprocess.call(remove, shell=True)
remove2='rm '+output+'/06_Final_Results/annotations/*.parsed.tsv'
subprocess.call(remove2, shell=True)
remove3='rm '+output+'/06_Final_Results/annotations/*.parsinator.csv'
subprocess.call(remove3, shell=True)
