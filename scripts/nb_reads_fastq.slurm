#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=500K
#SBATCH --time=00:30:00
#SBATCH --job-name=Quant_Analysis
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat/scripts/quant_analysis_%j.e


##########
# This script returns the total number of reads per FASTQ file
##########

# count the number of reads in the fastq files of parental cells and paraclonal cells for all replicates and store the result in the nb_reads_per_fastq.tsv file in the analysis/ directory  

for file in $@
do
	echo -e "File\t$file"
	gzip -dc $file | grep -c "@"
	echo ""
done >> nb_reads_per_fastq.tsv
