import sys,argparse,os,subprocess,re
import pandas as pd


##flags
parser = argparse.ArgumentParser(description='''Script to use PIGv in batch mode''')
parser.add_argument('-i', '--input', help='''Input folder of assemblies''', required=True)
parser.add_argument('-cov', '--coverm', help='''CoverM Folder (must have same filenames but .coverm rather than .fa)''', required=True)
parser.add_argument('-t', '--threads', help='''Number of threads to use''', required=True)
parser.add_argument('-annot', '--annotation', help='''Annotation mode mark true''', required=True)

args = parser.parse_args()
input_folder = args.input
coverm_folder = args.coverm
threads = args.threads
annotation = args.annotation

for file in os.listdir(input_folder):
    print('Running PIGv on '+file)
    assembly= os.path.join(input_folder, file)
    coverm_name=re.sub('.fa', '.coverm', file)
    coverm= os.path.join(coverm_folder, coverm_name)
    outname= re.sub('.contigs.fa','',file)
    outfolder= outname+'_PIGv'
    if annotation=='True':
        cmd1= 'python PIGv.py -i '+assembly+' -cov '+coverm+' -t '+threads+' -annot True -o '+ outfolder
    else:
        cmd1= 'python PIGv.py -i '+assembly+' -cov '+coverm+' -t '+threads+' -o '+ outfolder
    subprocess.call(cmd1, shell=True)
