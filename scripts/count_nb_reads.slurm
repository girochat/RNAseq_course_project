#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=500K
#SBATCH --time=00:30:00
#SBATCH --job-name=count_reads
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=quant_analysis_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# count the number of reads in each fastq file in argument
for file in $@
do
	file_name=$(basename $file)
	echo -e "File\t${file_name}"
	gzip -dc $file | grep -c "@"
	echo ""
done >> $HOME/analysis/results_nb_reads_fastq.tsv
