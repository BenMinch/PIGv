# PIGv with ViraScreen
![alt text](https://github.com/BenMinch/PIGv/blob/main/images/PIGv.png)
Pipeline for the identification of giant virus genomes from metagenomic datasets

# Description
An easy to install and use pipeline to bin, clean, annotate, classify, and quantify genomes of giant viruses from metagenomic datasets. This program also includes the software **ViraScreen** which can save you a lot of time by screening your raw fastq reads for signs of giant viruses before you start the pipeline. 

## Dependencies
It is recommended that you create a conda environment in which to install all of these packages to. 
1. Python packages (version 3.7.6)
- Pandas (version 1.3.5)
- natsort (version 8.3.1)
- Bio (version 1.5.8)
- numpy (version 1.21.5)
- matplotlib (version 3.5.3)
2. [HMMER3](https://github.com/EddyRivasLab/hmmer) version 3.3.2
3. [Prodigal-gv](https://github.com/apcamargo/prodigal-gv) version 2.11.0
4. Snakemake version 7.25
5. [Metabat2](https://bitbucket.org/berkeleylab/metabat) version 2.15
6. [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) version 1.0.1
7. [DIAMOND](https://github.com/bbuchfink/diamond/releases/tag/v2.0.4) version 2.0.4
7. [Megahit](https://github.com/voutcn/megahit) version 1.2.9
8. [CoverM](https://github.com/wwood/CoverM) version 0.4.0


## Installation

It's as easy as cloning this repository and making sure you have all the dependencies installed
`git clone https://github.com/BenMinch/PIGv`

### Downloading databases and setting up the program

Download the hmm database `wget -O hmm.tar.gz https://zenodo.org/record/4762520/files/hmm.tar.gz?download=1` and then `tar -xvzf hmm.tar.gz` 

Take all the files inside of this folder and copy them into the hmm folder provided inside the PIGv program from the github clone. 

Unpack the diamond program inside the viral_screening folder `gunzip diamond-linux64.tar.gz|tar -xvzf`
*Note: if you have a different computer model than linux64, you may need to install a different version of diamond and can do so [here](https://github.com/bbuchfink/diamond)*

Download the CheckV database and move it inside the resources folder `checkv download_database ./`

## ViraScreen

ViraScreen is a time-saver to scan your raw read files before you move forward with trimming, assembly and mapping in order to see if you have a good chance of getting a giant virus genome from the dataset. The program uses a blast search against a database of MCP proteins from 1500 known giant viruses and uses this data to predict the chance of a genome being inside your sample. This prediction is based on a model that was trained from over 50 uses of the program in conjunction with PIGv. 

### Usage

`python Viral_Screeninator.py -i input_folder`

1. -i: Input folder of raw fastq reads that you want to test. *Note: you may only want to test one of the two paired end reads for a speedier screening process* 

### Outputs

A single csv file is output that will contain columns for the filename, the number of MCP hits in that query file and the chance of getting a genome from this dataset (a number from 0-1 representing a proportion). 

# Using PIGv
![alt text](https://github.com/BenMinch/PIGv/blob/main/images/chart.png)
Once you have screened your reads and want to move forward with using PIGv, you'll find that it is quite easy. An example folder has some test inputs you can use to make sure the program is running correctly. It should give you 2 genomes. To run the examples you must first unzip the .fa files and run `cat korea.fa1 korea.fa2 > korea.fa` to combine the split fasta file (due to github filelimits).

**Different starting points**
1. Trimmed Reads: PIGv can do assembly and coverage mapping for you. All you need to input is trimmed forward and reverse reads (also works for single end reads as well). This mode will take significantly longer to run (about 2 hours more). 

*Minimal usage*
`python PIGv.py -1 read_1.fastq -2 read_2.fastq -o Output_directory -t 12 -type reads`

2. Assembled Contigs and CoverM file: If you already have assembled contigs and a coverm coverage file, you should run this mode as it is faster.

*Minimal usage*
`python PIGv.py -i assembled_contigs.fa -o Output_directory -t 12 -cov metagenome.coverm -annot True -type contigs`

**Inputs**
1. -i: assembled contigs (required if running type contigs)
2. -o: Output directory (can be whatever, it will create one if it doesn't exits) **Note: It cannot have "." or special characters. "_" is ok.
3. -t: threads. I usually use around 12
4. -cov: coverage file (required if running type contigs)
5. -annot: "True" if you want annotations for your genomes. This step does add around 30 minutes to completion time.
6. -1: fastq read fwd (required for type reads)
7. -2: fastq read rv (required for type reads)
8. -type: contig or reads mode (required)
9.-single: "True" if you have single end reads (you can input your read as the -1)
10. -contig: "True" if you want to do contig-level analysis for NCLDV Markersearch

**Outputs**
1. A folder of annotations (.annotated.csv) if you selected it
2. A folder of clean genomes containing your giant virus genomes
3. A file called gvclass_out.tab: this is the output from [gvclass](https://github.com/NeLLi-team/gvclass/) and it contains taxonomy information as well as other relevant genome statistics
4. A file called viral_genome_statistics.tsv: This contains information about different marker genes present in your genomes as well as the ViralRecall score. Read more about this [here](https://github.com/faylward/viralrecall)
5. If you want MCP hits you can find them in the folder labled Contig_markersearch in a table called mcp.table.tsv. The corresponding proteins for MCPs are in the mcp.faa file.
6. CheckV output file: file with information about completeness and quality of each genome.

## PIGv Batch

A script has been included to run PIGv on a plethora of files at the same time. The process is assentially the same but instead of an input assembly and coverage file, you can input a folder of assemblies and a separate folder of coverage files. This will create a separate output folder for each sample.
**NOTE:** this batch file requires that the name of your assemblies ends in ".contigs.fa" and that your coverage files end with ".contigs.coverm" with the basename being the same between the two files. This is how the program knows which ones to match.

PIGv Batch also works in reads mode and you can input a directory of trimmed fastq reads as the input. Just make sure the only differences between fwrd and reverse is "_1" and "_2".

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
