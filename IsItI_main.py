#!/usr/bin/python

#loading libraries
import sys
import argparse
import os
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com' 

def check_arg(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--single',
                        help='name of single-end fastq.gz file',
                        )
    parser.add_argument('-p1', '--pair1',
                        help='name of first paired-end fastq.gz file',
                        )
    parser.add_argument('-p2', '--pair2',
                        help='name of second paired-end fastq.gz file',
                        )
    parser.add_argument('-ts', '--txts',
                        help='name of text file with fastq.gz files',
                        )
    parser.add_argument('-r', '--reference',
                        help='name of reference fasta file',
                        )
    parser.add_argument('-tr', '--txtr',
                        help='name of text file with reference fasta files',
                        )
    parser.add_argument('-a', '--accession',
                        help='reference NCBI accession code',
                        )
    parser.add_argument('-ta', '--txta',
                        help='name of text file with reference NCBI accession codes',
                        )
    return parser.parse_args(args)

### GETTING INPUT FILES ###

# retrieve command line arguments
args = check_arg(sys.argv[1:])

#fastq files
if args.single != None: # check if single-end file was provided
    fastq_files = [args.single]
if args.pair1 != None and args.pair2: # check if paired-end file was provided
    fastq_files = [[args.pair1, args.pair2]]
if args.txts != None: # check if text file with fastq files was provided
    with open(args.txts) as input_fastq_file:
        lines = input_fastq_file.readlines()
        lines = [i.rstrip() for i in lines]
        lines = [i.split(' ') for i in lines]
        fastq_files = lines

#fasta files
if args.reference != None: # check if reference fasta file was provided
    fasta_files = [args.reference]
if args.txtr != None: # check if text file with reference fasta files was provided
    with open(args.txtr) as input_fasta_file:
        lines = input_fasta_file.readlines()
        lines = [i.rstrip() for i in lines]
        fasta_files = lines

#accession code
if args.accession != None: # check if reference accession code was provided
    handle = Entrez.efetch(db = 'nucleotide', id = args.accession, rettype = 'fasta', retmode = 'text') # retrieve sequence from NCBI
    record = handle.read()
    with open(args.accession+'.fasta', 'w') as reference_fasta: # save sequence into a file
        reference_fasta.write(record)
    fasta_files = [args.accession+'.fasta']
if args.txta != None: # check if text file with reference accession codes was provided
    with open(args.txta) as input_accession_file:
        fasta_files = []
        lines = input_accession_file.readlines()
        lines = [i.rstrip() for i in lines]
        for item in lines:
            handle = Entrez.efetch(db = 'nucleotide', id = item, rettype = 'fasta', retmode = 'text') # retrieve sequence from NCBI
            record = handle.read()
            with open(item+'.fasta', 'w') as reference_fasta: # save sequence into a file
                reference_fasta.write(record)
            fasta_files.append(item+'.fasta')

### RUNNING FASTQC ###

# checking if FastQC report folder exists. creates one if it doesnt. 
current_dir = os.getcwd()
if not os.path.isdir(current_dir+'/fastqc_reports'):
    os.makedirs(current_dir+'/fastqc_reports')
fastqc_reports_dir = current_dir+'/fastqc_reports'

# run FastQC for each fastq file, separately (even paired-ends)
for sequence in fastq_files:
    for i in range(len(sequence)):
        command_line = 'FastQC/fastqc -j ~/jdk-18/bin/java -o '+fastqc_reports_dir+' '+sequence[i]
        os.system(command_line)

### RUNNING TRIMMOMATIC ###

processed_fastq_files = []
for sequence in fastq_files:
    if len(sequence) == 1:
        sequence_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','')
        command_line = '~/jdk-18/bin/java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 '+sequence[0]+' '+sequence_name+'_out.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append(sequence_name+'_out.fq.gz')
    else:
        sequence0_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','')
        sequence1_name = sequence[1].replace('.fq.gz','').replace('.fastq.gz','')
        command_line = '~/jdk-18/bin/java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE '+sequence[0]+' '+sequence[1]+' '+sequence0_name+'_paired_out.fq.gz '+sequence0_name+'_unpaired_out.fq.gz '+sequence1_name+'_paired_out.fq.gz '+sequence1_name+'_unpaired_out.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append([sequence0_name+'_paired_out.fq.gz', sequence1_name+'_paired_out.fq.gz'])

### RUNNING FASTQC ON FILTERED SEQUENCES ###

# run FastQC for each fastq file, separately (even paired-ends)
for sequence in processed_fastq_files:
    for i in range(len(sequence)):
        command_line = 'FastQC/fastqc -j ~/jdk-18/bin/java -o '+fastqc_reports_dir+' '+sequence[i]
        os.system(command_line)