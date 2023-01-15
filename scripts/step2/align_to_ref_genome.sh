#!/usr/bin/env bash

module add UHTS/Aligner/hisat/2.2.1

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# add path of the genome index files to the HISAT2 environment variable 
HISAT2_INDEXES=/data/courses/rnaseq_course/lncRNAs/Project2/references/

# extract the generic name of the paired files for the output file
output_file=$(basename $1 | awk -F "_R" '{print $1}')

# save the result of the alignment from STDERR 
exec 2> $HOME/analysis/alignment/results_align_${output_file}.log

# align the paired-end reads to the human genome with HISAT2 
# option --dta : better compatibility with StringTie
# option --rna-strandness : first-strand stranded library
hisat2 -p 8 --rna-strandness RF --dta -x index_genome -1 $1 -2 $2 -S $HOME/data/${output_file}.sam
