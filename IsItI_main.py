#!/usr/bin/python

#loading libraries
import sys
import argparse
import os
from Bio import Entrez
from Bio import SeqIO
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

log_output = open('IsItI_log.txt', 'w') # create and open log file

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
log_output.write('Initial sample file quality assessment successfully finalized. Results can be found in the "fastqc_reports" folder\n') # write into log file

### RUNNING TRIMMOMATIC ###

processed_fastq_files = []
for sequence in fastq_files: # for every fastq item file in the list
    if len(sequence) == 1: # check if it is single-end
        sequence_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        command_line = 'jdk-18/bin/java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 '+sequence[0]+' '+sequence_name+'_out.fq.gz ILLUMINACLIP:'+current_dir+'/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append([sequence_name+'_out.fq.gz']) # append output file name to list
        log_output.write(sequence_name+" successfully QC'd.\n") # write into log file
    else: # else, it is paired-end
        sequence0_name = sequence[0].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        sequence1_name = sequence[1].replace('.fq.gz','').replace('.fastq.gz','') # get file name without the extension
        command_line = 'jdk-18/bin/java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 '+sequence[0]+' '+sequence[1]+' '+sequence0_name+'_out.fq.gz '+sequence0_name+'_unpaired_out.fq.gz '+sequence1_name+'_out.fq.gz '+sequence1_name+'_unpaired_out.fq.gz ILLUMINACLIP:'+current_dir+'/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        os.system(command_line)
        processed_fastq_files.append([sequence0_name+'_out.fq.gz', sequence1_name+'_out.fq.gz']) # append output file name to list
        log_output.write(sequence0_name+" and "+sequence1_name+" successfully QC'd.\n") # write into log file

### RUNNING FASTQC ON FILTERED SEQUENCES ###

# run FastQC for each fastq file, separately (even paired-ends)
for sequence in processed_fastq_files:
    for i in range(len(sequence)):
        command_line = 'FastQC/fastqc -j jdk-18/bin/java -o '+fastqc_reports_dir+' '+sequence[i]
        os.system(command_line)
log_output.write('Sample file quality assessment after QC successfully finalized. Results can be found in the "fastqc_reports" folder\n') # write into log file

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
sam_files = [] # list to store sam files names
sample_to_reference = {} # create a dictionary to keep track of sample-reference pairs
for sequence in processed_fastq_files: # for each sample fastq file (QC'd)
    for reference in fasta_files: #  for each reference genome
        seq_file_name = sequence[0].replace('_out.fq.gz','') # remove extension to get file name
        ref_file_name = reference.replace('.fasta','').replace('.fa','') # remove extension to get file name

        if seq_file_name not in sample_to_reference: # if sample file is not in the sample_to_reference dict
            sample_to_reference[seq_file_name]=[ref_file_name] # add the new sample-reference pair
        else: # if it already exists
            sample_to_reference[seq_file_name].append(ref_file_name) # just add a new reference to the values

        if len(sequence) == 1: # check if fastq file is single-end
            command_line = 'bowtie2-2.4.5-linux-x86_64/bowtie2 -x '+bowtie_files_dir+'/'+ref_file_name+' -U '+sequence[0]+' -S '+seq_file_name+'_mappedto_'+ref_file_name+'.sam'
            sam_files.append(seq_file_name+'_mappedto_'+ref_file_name+'.sam') # append output SAM file name to list
            os.system(command_line) # runs bowtie2
        else: # else, it is paired-end
            command_line = 'bowtie2-2.4.5-linux-x86_64/bowtie2 -x '+bowtie_files_dir+'/'+ref_file_name+' -1 '+sequence[0]+' -2 '+sequence[1]+' -S '+seq_file_name+'_mappedto_'+ref_file_name+'.sam'
            sam_files.append(seq_file_name+'_mappedto_'+ref_file_name+'.sam') # append output SAM file name to list
            os.system(command_line) # runs bowtie2
        log_output.write(seq_file_name+' successfully mapped to '+ref_file_name+'\n') # write into log file

### Running SAMTOOLS to DETERMINE Sequence Coverage of Reference Genome ###

# preparing input files for mpileup
bam_files = [] # create a new list that will contain converted BAM files
for sam_file in sam_files: # iterate through each SAM file and convert to BAM file
    bam_file = sam_file.replace('.sam','.bam') # replace extension
    sam_to_bam = './samtools-1.9/samtools view -Sb '+sam_file+' > '+bam_file # command to convert sam file to bam file 
    sort_bam = './samtools-1.9/samtools sort '+bam_file+' -o '+bam_file.replace('.bam','.sort.bam') # command to sort bam file 
    index_bam = './samtools-1.9/samtools index '+bam_file.replace('.bam','.sort.bam') # command to index bam file 
    bam_files.append(bam_file.replace('.bam','.sort.bam')) # append indexed & sorted bam file to list
    os.system(sam_to_bam) # os call to run sam_to_bam command
    os.system(sort_bam) # os call to run sort_bam command
    os.system(index_bam) # os call to run index_bam command

# run mpileup
for bam_file in bam_files: # for each bam file stored in bam_files
    reference_name = bam_file[bam_file.find('mappedto_')+len('mappedto_'):bam_file.find('.sort.bam')]+'.fasta' # get the correct reference file name
    if not os.path.exists(reference_name): # in case the reference fasta file does not end with .fasta
        reference_name = bam_file[bam_file.find('mappedto_')+len('mappedto_'):bam_file.find('.sort.bam')]+'.fa' # change it to .fa
    mpileup_command = './samtools-1.9/samtools mpileup -f '+reference_name+' '+bam_file+' > '+bam_file.replace('.sort.bam','.mpileup') # create mpileup command 
    os.system(mpileup_command) # os call to run mpileup command
    
# run flagstat
for bam_file in bam_files: # for each bam file stored in bam_files
    flagstat_command = './samtools-1.9/samtools flagstat '+bam_file+' > '+bam_file.replace('.sort.bam','.flagstat') # create the samtools flagstat command
    os.system(flagstat_command) #write command out to os.system

### SUMMARIZING SAMTOOLS RESULTS ###
 
# we need to know the size of each reference genome to compute average coverage
reference_sizes = {} # empty dictionary to store reference id and length 
for item in fasta_files: # for each reference file
    record = list(SeqIO.parse(item, 'fasta')) # read fasta file
    for sequence in record: # read iterator
        reference_sizes[sequence.id[:sequence.id.find('.')]]=len(str(sequence.seq)) # add to dictionary reference id as key, sequence length as value

# initialize final output file
summarized_output = open('IsItI_results.csv', 'w') # create and open file
summarized_output.write('Sample,Reference,Reads in Sample,Reads Mapped to Reference,Average Coverage (reads/nt)\n') # write file header

for sample in sample_to_reference: # for each sample
    for reference in sample_to_reference[sample]: # for each reference mapped to
        with open(sample+'_mappedto_'+reference+'.flagstat') as flagstat_file: # open corresponding flagstat file
            for line in flagstat_file:
                if line.find('in total') != -1: # find line with total number of reads count
                    n_reads_sample = int(line[:line.find(' + ')]) # save that number 
                if line.find('mapped') != -1 and line.find('%') != -1: # find line with number of reads mapped to reference
                    n_reads_map = int(line[:line.find(' + ')]) # save that number

        if n_reads_map == 0: # if the number of reads mapped is zero, no need to open mpileup file as average coverage will be zero
            summarized_output.write(sample+','+reference+','+str(n_reads_sample)+','+str(n_reads_map)+',0\n') # write results into file
        else: # if at least one read mapped 
            with open(sample+'_mappedto_'+reference+'.mpileup') as mpileup_file: # open corresponding mpileup file
                reads_mapped = 0 # set initial value to zero
                for line in mpileup_file: # for each line in file
                    line = line.split() # split file into a list
                    reads_mapped += int(line[3]) # get the element on the fourth column, which corresponds to the number of reads mapped to a specific position
            avg_cov = reads_mapped/reference_sizes[reference] # compute average coverage
            summarized_output.write(sample+','+reference+','+str(n_reads_sample)+','+str(n_reads_map)+','+str(avg_cov)+'\n') # write results into file
        
        log_output.write('Mapping results for '+sample+' and '+reference+' have been summarized into "IsItI_results.csv"\n') # write into log file

# closing all files
log_output.close()
summarized_output.close()