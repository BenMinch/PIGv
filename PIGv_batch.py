import sys,argparse,os,subprocess,re
import pandas as pd


##flags
parser = argparse.ArgumentParser(description='''Script to use PIGv in batch mode''')
parser.add_argument('-i', '--input', help='''Input folder of assemblies''', required=False)
parser.add_argument('-ref', '--reference', help='''Reference CSV of assemblies and fastq names [for contig mode]''', required=False)
parser.add_argument('-reads', '--reads', help='''Reads folder ''', required=False)
parser.add_argument('-t', '--threads', help='''Number of threads to use''', required=True)
parser.add_argument('-annot', '--annotation', help='''Annotation mode mark true''', required=True)
parser.add_argument('-type', '--type',help='Type of input: contigs or reads', choices=['contigs', 'reads'], default='contigs')
parser.add_argument('-n', '--name', help='Name of the run', required=True)
args = parser.parse_args()
input_folder = args.input
read_folder = args.reads
reference= args.reference
threads = args.threads
annotation = args.annotation
name=args.name
type = args.type


if type == 'contigs':
    print('Preparing CoverM files')
    coverage= 'python scripts/Mapping_Pipeline_metagenome_binning.py '+read_folder+' '+'coverage'+' '+input_folder+' '+reference
    subprocess.call(coverage, shell=True)
    for file in os.listdir(input_folder):
        print('Running PIGv on '+file)
        assembly= os.path.join(input_folder, file)
        coverm_name=re.sub('.fa', '.coverm', file)
        coverm= os.path.join('coverage', coverm_name)
        outname= re.sub('_contigs.fa','',file)
        outfolder= outname+'_PIGv'
        if annotation=='True':
            cmd1= 'python PIGv.py -i '+assembly+' -cov '+coverm+' -t '+threads+' -annot True -type contigs -o '+ name+'_'+outfolder+' -n '+name+'_'+outname
        else:
            cmd1= 'python PIGv.py -i '+assembly+' -cov '+coverm+' -t '+threads+' -type contigs -o '+ name+'_'+outfolder + ' -n '+name+'_'+outname
        subprocess.call(cmd1, shell=True)

else:
    for file in os.listdir(read_folder):
        print('Running PIGv on '+file)
        R1= os.path.join(read_folder, file)
        R2= re.sub('_1', '_2', R1)
        if R1==R2:
            pass
        else:
            outname= re.sub('.fq.gz','',file)
            outfolder= outname+'_PIGv'
            if annotation=='True':
                cmd1= 'python PIGv.py -1 '+R1+' -2 '+R2+' -t '+threads+' -annot True -type reads -o '+ name+'_'+outfolder+' -n '+name+'_'+outname
            else:
                cmd1= 'python PIGv.py -1 '+R1+' -2 '+R2+' -t '+threads+' -type reads -o '+ name+'_'+outfolder + ' -n '+name+'_'+outname

            subprocess.call(cmd1, shell=True)
  
            
