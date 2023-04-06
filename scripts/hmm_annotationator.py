import pandas as pd
import numpy as np
import os, argparse, sys, re, subprocess
#### HMM annotationinator
# This script takes a parsed HMM output file and a reference annotation and combines them

# define flags
parser = argparse.ArgumentParser(description='A script to combine HMM output and reference annotation')
parser.add_argument('-i', '--input', help='Input HMM output file', required=True)
parser.add_argument('-r', '--reference', help='Reference annotation file', required=True)
parser.add_argument('-o', '--output', help='Output file name', required=True)

args = parser.parse_args()

# set variables
input_file = args.input
reference_file = args.reference
output_file = args.output

# read in files
input_df = pd.read_csv(input_file)
reference_df = pd.read_csv(reference_file,sep='\t')

#add columns to input_df
input_df['NCVOG_descs']= input_df['id'].map(reference_df.set_index('GVOG')['NCVOG_descs'])
input_df['consensus_NOGs']= input_df['id'].map(reference_df.set_index('GVOG')['consensus_NOGs'])
input_df['nog_descs']= input_df['id'].map(reference_df.set_index('GVOG')['nog_descs'])
input_df['nog_categories']= input_df['id'].map(reference_df.set_index('GVOG')['nog_categories'])
input_df['pfam_descs']= input_df['id'].map(reference_df.set_index('GVOG')['pfam_descs'])

# write output
input_df.to_csv(output_file, index=False)

