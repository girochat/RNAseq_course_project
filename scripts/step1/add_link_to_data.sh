#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# path to the reference files for the project
REF=/data/courses/rnaseq_course/lncRNAs/fastq

# add symbolic link to the fastq files of the paraclonal cells in data/ directory
for i in 2 4 7
do
	ln -s $REF/3_${i}_L3_R1_*.fastq.gz $HOME/data/3_${i}_R1.fastq.gz
	ln -s $REF/3_${i}_L3_R2_*.fastq.gz $HOME/data/3_${i}_R2.fastq.gz
done

# add symbolic link to the fastq files of the parental cells in data/ directory
for i in 1 2 3
do
	ln -s $REF/P${i}_L3_R1_*.fastq.gz $HOME/data/P${i}_R1.fastq.gz
	ln -s $REF/P${i}_L3_R2_*.fastq.gz $HOME/data/P${i}_R2.fastq.gz
done
