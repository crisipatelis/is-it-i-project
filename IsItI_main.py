#!/usr/bin/python

#loading libraries
import sys
import argparse
import os
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com' 

def check_arg(args=None): # define possible parameters that can be used
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--single',
                        help='name of single-end fastq.gz file')
    parser.add_argument('-p1', '--pair1',
                        help='name of first paired-end fastq.gz file')
    parser.add_argument('-p2', '--pair2',
                        help='name of second paired-end fastq.gz file')
    parser.add_argument('-ts', '--txts',
                        help='name of text file with fastq.gz files')
    parser.add_argument('-r', '--reference',
                        help='name of reference fasta file')
    parser.add_argument('-tr', '--txtr',
                        help='name of text file with reference fasta files')
    parser.add_argument('-a', '--accession',
                        help='reference NCBI accession code')
    parser.add_argument('-ta', '--txta',
                        help='name of text file with reference NCBI accession codes')
    return parser.parse_args(args)

### GETTING INPUT FILES ###

# retrieve command line arguments
args = check_arg(sys.argv[1:])

# checking if all required inputs were given
if not (args.single != None or args.pair1 != None and args.pair2 != None or args.txts != None): # check sample input
    print('Sample input not informed. At least one FASTQ input file is required to run this tool.')
    exit()
if not (args.reference != None or args.txtr != None or args.accession != None or args.txta != None): # check reference input
    print('Reference input not informed. At least one FASTA file or Accession Code is required to run this tool.')
    exit()

#fastq files
if args.single != None: # check if single-end file was provided
    fastq_files = [args.single] # append file name to list
elif args.pair1 != None and args.pair2 != None: # check if paired-end file was provided
    fastq_files = [[args.pair1, args.pair2]]  # append list of file names to list (yes, a list inside of a list)
elif args.txts != None: # check if text file with fastq files was provided
    with open(args.txts) as input_fastq_file:
        lines = input_fastq_file.readlines() # read text file into a list 
        lines = [i.rstrip() for i in lines] # remove newline character
        lines = [i.split(' ') for i in lines] # if paired end files are in the text file, splits item in list by space 
        fastq_files = lines 

#fasta files
if args.reference != None: # check if reference fasta file was provided
    fasta_files = [args.reference] # append file name to list
elif args.txtr != None: # check if text file with reference fasta files was provided
    with open(args.txtr) as input_fasta_file:
        lines = input_fasta_file.readlines() # read text file into a list 
        lines = [i.rstrip() for i in lines] # remove newline character
        fasta_files = lines

#accession code
if args.accession != None: # check if reference accession code was provided
    handle = Entrez.efetch(db = 'nucleotide', id = args.accession, rettype = 'fasta', retmode = 'text') # retrieve sequence from NCBI
    record = handle.read() # read it
    with open(args.accession+'.fasta', 'w') as reference_fasta: # save sequence into a file
        reference_fasta.write(record)
    fasta_files = [args.accession+'.fasta'] # append file name to list
elif args.txta != None: # check if text file with reference accession codes was provided
    with open(args.txta) as input_accession_file:
        fasta_files = []
        lines = input_accession_file.readlines() # read text file into a list of accession numbers
        lines = [i.rstrip() for i in lines] # remove newline character
        for item in lines:
            handle = Entrez.efetch(db = 'nucleotide', id = item, rettype = 'fasta', retmode = 'text') # retrieve sequence from NCBI
            record = handle.read()
            with open(item+'.fasta', 'w') as reference_fasta: # save sequence into a file
                reference_fasta.write(record)
            fasta_files.append(item+'.fasta') # append file name to list

### RUNNING FASTQC ###

# checking if FastQC report folder exists. creates one if it doesn't. 
current_dir = os.getcwd()
if not os.path.isdir(current_dir+'/fastqc_reports'):
    os.makedirs(current_dir+'/fastqc_reports')
fastqc_reports_dir = current_dir+'/fastqc_reports'

# run FastQC for each fastq file, separately (even paired-ends)
for sequence in fastq_files:
    for i in range(len(sequence)):
        command_line = 'FastQC/fastqc -j jdk-18/bin/java -o '+fastqc_reports_dir+' '+sequence[i]
        os.system(command_line)

### RUNNING TRIMMOMATIC ###

processed_fastq_files = []
for sequence in fastq_files: # for every fastq item file in the list
    if len(sequence) == 1: # check if it is single-end
        sequence_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        command_line = 'jdk-18/bin/java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 '+sequence[0]+' '+sequence_name+'_out.fq.gz ILLUMINACLIP:'+current_dir+'/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append([sequence_name+'_out.fq.gz']) # append output file name to list
    else: # else, it is paired-end
        sequence0_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        sequence1_name = sequence[1].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        command_line = 'jdk-18/bin/java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 '+sequence[0]+' '+sequence[1]+' '+sequence0_name+'_out.fq.gz '+sequence0_name+'_unpaired_out.fq.gz '+sequence1_name+'_out.fq.gz '+sequence1_name+'_unpaired_out.fq.gz ILLUMINACLIP:'+current_dir+'/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append([sequence0_name+'_out.fq.gz', sequence1_name+'_out.fq.gz']) # append output file name to list

### RUNNING FASTQC ON FILTERED SEQUENCES ###

# run FastQC for each fastq file, separately (even paired-ends)
for sequence in processed_fastq_files:
    for i in range(len(sequence)):
        command_line = 'FastQC/fastqc -j jdk-18/bin/java -o '+fastqc_reports_dir+' '+sequence[i]
        os.system(command_line)

### RUNNING BOWTIE2 TO MAP READS TO REFERENCE GENOME ###

# checking if bowtie folder exists. creates one if it doesn't. 
if not os.path.isdir(current_dir+'/bowtie_files'):
    os.makedirs(current_dir+'/bowtie_files')
bowtie_files_dir = current_dir+'/bowtie_files'

# create a bowtie2 index to map to for each reference fasta file
for reference_file in fasta_files: # for each fasta file in the list
    file_name = reference_file.replace('.fasta','').replace('.fa','') # remove extension to get file name
    command_line = 'bowtie2-2.4.5-linux-x86_64/bowtie2-build '+reference_file+' '+bowtie_files_dir+'/'+file_name
    os.system(command_line) # runs bowtie2-build

# mapping fastq files to reference genomes 
sam_files = [] 
sequence_names = []
for sequence in processed_fastq_files: # for each sample fastq file (QC'd)
    for reference in fasta_files: #  for each reference genome
        seq_file_name = sequence[0].replace('_out.fq.gz','') # remove extension to get file name
        sequence_names.append(seq_file_name) #names for the samtools output
        ref_file_name = reference.replace('.fasta','').replace('.fa','') # remove extension to get file name
        if len(sequence) == 1: # check if fastq file is single-end
            command_line = 'bowtie2-2.4.5-linux-x86_64/bowtie2 -x '+bowtie_files_dir+'/'+ref_file_name+' -U '+sequence[0]+' -S '+seq_file_name+'_mappedto_'+ref_file_name+'.sam'
            sam_files.append(seq_file_name+'_mappedto_'+ref_file_name+'.sam') # append output SAM file name to list
            os.system(command_line) # runs bowtie2
        else: # else, it is paired-end
            command_line = 'bowtie2-2.4.5-linux-x86_64/bowtie2 -x '+bowtie_files_dir+'/'+ref_file_name+' -1 '+sequence[0]+' -2 '+sequence[1]+' -S '+seq_file_name+'_mappedto_'+ref_file_name+'.sam'
            sam_files.append(seq_file_name+'_mappedto_'+ref_file_name+'.sam') # append output SAM file name to list
            os.system(command_line) # runs bowtie2

### Running SAMTOOLS to DETERMINE Sequence Coverage of Reference Genome ###

for i in range(len(sam_files)): #run through all SAM files
    sam_flagstat_command = ‘samtools flagstat ‘ + sam_files[i] + ‘ -o current_dir sam_’ + sequence_names[i] + ‘_statsout.tsv’ #create the samtools flagstat command
    os.system(sam_flagstat_command) #write command out to os.system
