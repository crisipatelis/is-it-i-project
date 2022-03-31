<!-----

Yay, no errors, warnings, or alerts!

Conversion time: 0.459 seconds.


Using this Markdown file:

1. Paste this output into your source file.
2. See the notes and action items below regarding this conversion run.
3. Check the rendered output (headings, lists, code blocks, tables) for proper
   formatting and use a linkchecker before you publish this page.

Conversion notes:

* Docs to Markdown version 1.0β33
* Thu Mar 31 2022 09:29:36 GMT-0700 (PDT)
* Source doc: COMP 483 Design Doc
----->


**A Mapping-Based Approach to Identify Species of Interest**

**Daniel Araujo, Michael Lemenager, Crisi Patelis**

**Overview** 

Often scientists are interested in organisms that do not properly grow in lab conditions and/or by themselves, and therefore are more difficult to have their genomes sequenced and analyzed. To solve this issue, metagenomics approaches, such as metabarcoding and shotgun, have been widely implemented to sequence all (or most of) DNA contained within a sample. 

A crucial step after sequencing DNA in a sample is to determine which species sequencing reads belong to. This is particularly difficult with short reads sequencing, which is the most commonly implemented method, as short reads do not possess enough information by themselves to identify organisms in a sample. Thus, a common step is to assemble reads into larger contigs, which are more informative. However, assembling complete genomes of all organisms contained within a sample is not an easy task, as genome assembly can be influenced by sequencing coverage, presence of closely-related species in a sample, repetitive regions on the DNA, and other intrinsic features of the sample. Therefore, it is necessary to have different approaches not based on genome assemblies to identify species of interest in a sample. One alternative way would be to map short sequencing reads to reference genomes or gene markers. 

**Context** 

          When identifying a species of interest within a given sample, there are two particular challenges that arise. The first challenge is that genomes cannot be assembled for each organism in order to identify them. With the number of organisms that may be contained within a particular sample, it is both inefficient and infeasible to sequence the genomes of each and compare them to a reference genome. The second problem is that closely related species need to be differentiated, which can be challenging when they have similar genomic compositions. 

	Our pipeline is necessary because it aims to combat these challenges. The project will map reads to a given reference genome from a species of interest or to gene markers. In doing so, shorter reads are compared to the larger reference genome to analyze how much coverage matches. The user will be able to identify if their species of interest is present if enough reads are mapped against the genome. This solves a way to both compare sequence reads to a species of interest’s genome (without assembling all other genomes) and distinguish it from other closely related species based on a required amount of mapping coverage. 

	

**Goals & Non-Goals**

The goals for this project are:



* Develop a pipeline in Python to identify whether species of interest are present in NGS fastq files.
* For a more user-friendly implementation, this pipeline along with all its dependencies will be contained within a Docker image container. 

The non-goals for this project are:



* This project does not aim to build new assembly or mapping softwares. 

**Milestones**



* Monday, March 21st:
    * Design document should be finished.
    * Collection of test data to use during pipeline development. 
* Monday, March 28th:
    * Selection of fastq quality control software, and metrics to filter by. 
    * Pipeline should provide library quality reports, before and after QC, for each fastq file provided by the user. 
    * Initial commit to the GitHub repository.
* Monday, April 4th:
    * Mapping tool selected and tested to ensure that it is running correctly and giving the expected output.
    * Ensure that upstream processes are fully integrated to the mapping step. Report any errors found and solve them. 
    * Update code in the GitHub repository, and write initial user documentation (README file).  
* Monday, April 11th: 
    * Ensure that the pipeline is reporting proper mapping metrics, such as percentage of reads mapped to the reference genomes and overall coverage.
    * Familiarization with Docker image container building. 
    * Ensure that all processes so far are fully integrated. Report any errors found and solve them. 
    * Update code in the GitHub repository, and user documentation. 
* Monday, April 18th:
    * Ensure that all processes so far are fully integrated. Report any errors found and solve them. 
    * Start building Docker image container. Troubleshoot any problems found. 
    * Update code in the GitHub repository, and user documentation. 
* Monday, April 25th: 
    * The pipeline should be tested in different scenarios to confirm that it is running as expected, without any errors. Any unexpected outcome will be reported and solved. 
    * Metrics, such as total runtime of each scenario, will be reported for proper documentation. 
    * Final code and user documentation should be uploaded to the GitHub repository.
    * Docker image container should be finished.
* Monday, May 2nd: 
    * Project complete.

**Proposed Solution**

The following proposed solution is not final and can be modified if necessary. 



1. User inputs:
    1. The user will provide fastq files, which may or may not contain sequences derived from the species of interest.
        1. We expect that often users will provide raw fastq files. Thus, they should undergo QC processing, such as removal of adapter sequences and filtering low quality reads (suggestion of tool to use: Trimmomatic, BBDuk). 
        2. We may want to assess quality parameters, such as number of reads and overall read quality before and after QC, and return this to the user (using FastQC).
        3. The pipeline we design should be able to handle taking multiple fastq files as input. For this, the user will either supply the name of a single fastq file, or a text file containing the names of multiple fastq files - one file name per line. If the user provides more than one fastq file, they will be processed one at a time. If there are paired samples, file names should be on the same line in the text file, separated by a space character. 
    2. The user will also provide fasta files or accession codes corresponding to the species of interest.
        4. There will be no processing of sequences in the fasta files; they will be used as is.
        5. The pipeline should be able to handle multi-fasta files as well, as long as all sequences within the multi-fasta file belong to the same organism (multiple contigs).
        6. The pipeline we design should be able to handle taking multiple fasta files as input. For this, the user will either supply the name of a single fasta file, or a text file containing the names of multiple fasta files - one file name per line. If the user provides more than one fasta file, they will be processed one at a time.  
        7. If accession codes are provided instead of fasta files, the handling process will be similar to the described on 1.b.iii - the user can either supply a single accession code, or a text file containing multiple accession codes, one per line. Accession codes supplied should match NCBI’s database. 
2. Assessing presence of species of interest:
    3. As previously mentioned, we will not be relying on assembly-based methods to identify species of interest. Rather, we will be using a mapping approach. 
        8. After the initial QC, processed fastq reads will be mapped to the genomes of interest provided by the user. For this, we will look for in the literature benchmarking studies of mapping tools and identify which one is the most commonly recommended based on less memory usage, faster x, and higher accuracy. Suggestions include Bowtie2, BWA, and SOAP2. 
        9. If fastq files provided are paired ends, they will be aligned first, and then mapped to the reference genomes.
        10. We will assess important statistics such as number of reads mapped to the genome (in percentage) and mapping coverage of the genome using SAMtools, and provide those metrics to the user. The two aforementioned metrics are essential because by analyzing them, the user will know whether the species of interest are in the fastq files provided. 
        11. If multiple fastq files were provided, and if multiple reference fasta files/accession codes were provided as well, each fastq-fasta pair will be analyzed individually, but at the end, all results (mentioned in 2.a.ii) will be summarized in a single output file to facilitate user interpretation. 
3. Portability of our tool:
    4. Our pipeline, alongside all its dependencies, will be contained within a Docker image container. This way, users can use it regardless of operating system. This also facilitates the usage of our tool by those who are not familiar with bioinformatics analyzes. 

A graphical representation of our proposed solution is included below:

