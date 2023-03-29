import os
import sys
import subprocess
import re
import shlex

input_folder = sys.argv[1]
output_folder= sys.argv[2]
extension= sys.argv[3]
for i in os.listdir(input_folder):
	folder= os.path.join(input_folder,i)
	cmd1= 'cp '+ folder + '/*'+ extension+' '+ output_folder
	subprocess.call(cmd1, shell=True)

