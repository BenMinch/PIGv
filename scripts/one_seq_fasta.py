#Combine contigs into one sequence
import os,sys,argparse

def combine_contigs(file_path):
    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Get the file name without the file extension
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        # Initialize a variable to store the combined sequence
        combined_seq = ""
        # Iterate through the file lines
        for line in file:
            # If the line starts with '>', it's a header
            if line.startswith(">"):
                # Skip the line
                continue
            # Otherwise, it's a sequence line
            else:
                # Add the line to the combined sequence
                combined_seq += line.strip()
    # Write the combined sequence to a new file
    with open(file_name+'_combined.fasta', 'w') as new_file:
        new_file.write(">" + file_name + "\n")
        new_file.write(combined_seq)


file=sys.argv[1]
combine_contigs(file)