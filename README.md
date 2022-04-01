# Is it I? - A Mapping-Based Approach to Identify Species of Interest in FASTQ Files

## Introduction 

Often we are interested in organisms that do not properly grow in lab conditions and/or by themselves, and therefore are more difficult to have their genomes sequenced and analyzed. To solve this issue, metagenomics approaches, such as metabarcoding and shotgun, have been widely implemented to sequence all (or most of) DNA contained within a sample. A crucial step after sequencing DNA in a sample is to determine which species sequencing reads belong to, and a common step is to assemble reads into larger contigs. However, assembling complete genomes of all organisms contained within a sample is not an easy task, as genome assembly can be influenced by sequencing coverage, presence of closely-related species in a sample, repetitive regions on the DNA, and other intrinsic features of the sample. Therefore, it is necessary to have different approaches not based on genome assemblies to identify species of interest in a sample.

Our `Is it I?` pipeline avoids those challanges because it follows a different approach. Instead of assembly genomes, it maps reads to a given reference genome from a species of interest or to gene markers. In doing so, shorter reads are compared to the larger reference genome to analyze how much coverage matches. The user will be able to identify if their species of interest is present if enough reads are mapped against the genome. This solves a way to both compare sequence reads to a species of interest’s genome (without assembling all other genomes) and distinguish it from other closely related species based on a mapping coverage threshold.

## Implementation

Our `Is it I?` pipeline is written entirely in Python3 to be used in Linux OS, and uses well-maintained bioinformatics tools. Furthermore, our pipeline will be available at Docker Hub Repository, which means that users can use our tool regardless of OS, withouting having to worry about installing all dependencies. 

## Instructions

(To be added.)

## Credits

Daniel Araujo, Michael Lemenager, Crisi Patelis
