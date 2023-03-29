# PIGv
Pipeline for the identification of giant virus genomes from metagenomic datasets

# Description
An easy to install and use pipeline to bin, clean, annotate, classify, and quantify genomes of giant viruses from metagenomic datasets. This program also includes the software **ViraScreen** which can save you a lot of time by screening your raw fastq reads for signs of giant viruses before you start the pipeline. 

## Dependencies

1. Python packages
- Pandas
- Bio
- numpy
- matplotlib
2. [HMMER3](https://github.com/EddyRivasLab/hmmer)
3. [Prodigal-gv](https://github.com/apcamargo/prodigal-gv)
4. Snakemake
5. [Metabat2](https://bitbucket.org/berkeleylab/metabat)

## Installation

It's as easy as cloning this repository
`git clone https://github.com/BenMinch/PIGv`

## Preprocessing raw fastq

The inputs for the PIGv program are (1) a file of assembled contigs from a metagenomic dataset, (2) a coverage file in metabat2 format. This is how I go about getting these files from raw reads.

1. Trim the reads with trimgalore or another equivalent
2. Assemble the reads using megahit or another assembly program.
3. Map the trimmed reads back onto the assembly using CoverM contig with the method "metabat". You can read how to do this [Here](https://github.com/wwood/CoverM)

## ViraScreen

ViraScreen is a time-saver to scan your raw read files before you move forward with trimming, assembly and mapping in order to see if you have a good chance of getting a giant virus genome from the dataset. The program uses a blast search against a database of MCP proteins from 1500 known giant viruses and uses this data to predict the chance of a genome being inside your sample. This prediction is based on a model that was trained from over 50 uses of the program in conjunction with PIGv. 

### Usage

`python Viral_Screeninator.py -i input_folder`

1. -i: Input folder of raw fastq reads that you want to test. *Note: you may only want to test one of the two paired end reads for a speedier screening process* 

### Outputs

A single csv file is output that will contain columns for the filename, the number of MCP hits in that query file and the chance of getting a genome from this dataset (a number from 0-1 representing a proportion). 

# Using PIGv
Once you have screened your reads and want to move forward with using PIGv, you'll find that it is quite easy. An example folder has some test inputs you can use to make sure the program is running correctly. It should give you 2 genomes. 

`python PIGv.py -i assembled_contigs.fa -o Output_directory -t 12 -cov metagenome.coverm -annot True`

**Inputs**
1. -i: assembled contigs (required)
2. -o: Output directory (can be whatever, it will create one if it doesn't exits) **Note: It cannot have "." or special characters. "_" is ok.
3. -t: threads. I usually use around 12
4. -cov: coverage file (required)
5. -annot: "True" if you want annotations for your genomes. This step does add around 30 minutes to completion time.

**Outputs**
1. A folder of annotations (.annotated.csv) if you selected it
2. A folder of clean genomes containing your giant virus genomes
3. A file called gvclass_out.tab: this is the output from [gvclass](https://github.com/NeLLi-team/gvclass/) and it contains taxonomy information as well as other relevant genome statistics
4. A file called viral_genome_statistics.tsv: This contains information about different marker genes present in your genomes as well as the ViralRecall score. Read more about this [here](https://github.com/faylward/viralrecall)

## PIGv Batch

A script has been included to run PIGv on a plethora of files at the same time. The process is assentially the same but instead of an input assembly and coverage file, you can input a folder of assemblies and a separate folder of coverage files. This will create a separate output folder for each sample.
**NOTE:** this batch file requires that the name of your assemblies ends in ".contigs.fa" and that your coverage files end with ".contigs.coverm" with the basename being the same between the two files. This is how the program knows which ones to match. 

## References and Resources

This wrapper tool is built on the foundation of many other great tools.
- [Metabat2](https://bitbucket.org/berkeleylab/metabat)
- [ViralRecall](https://github.com/faylward/viralrecall)
- [NCLDV_markersearch](https://github.com/faylward/ncldv_markersearch)
- [Gvclass](https://github.com/NeLLi-team/gvclass/)

Screening of giant viruses was adapted from [Aylward et al. 2021](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001430).

# Copywright

PIGv Copyright (C) 2023 Benjamin Minch 

This program is free software: you can redistribute it and/or modify it under the terms of the MIT License. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details. 
