import os
import sys
import subprocess
import re
import shlex

input_folder = sys.argv[1]
outfolder = sys.argv[2]
i=input_folder
name = re.sub("_genomic.fna", "", i)
fasta = i
faa = os.path.join(outfolder,"prodigal.faa")
gff = os.path.join(outfolder, name+".gff")
genes = os.path.join(outfolder, name+".genes.fna")

cmd = "scripts/prodigal-gv -i "+ fasta +" -a "+ faa +' -p meta'
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))


