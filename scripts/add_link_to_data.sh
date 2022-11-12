#!/usr/bin/env bash

#################
# This script adds symbolic link to the FASTQ files of the cell subpopulations under analysis
################

# add symbolic link to the fastq files of the meroclonal cells in data/ directory
for i in 2 4 7
do
ln -s /data/courses/rnaseq_course/lncRNAs/fastq/3_${i}_L3_R1_*.fastq.gz .
done

# add symbolic link to the fastq files of the parental cells in data/ directory
for i in 1 2 3
do
ln -s /data/courses/rnaseq_course/lncRNAs/fastq/P${i}_L3_R1_*.fastq.gz .
done
