#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=250M
#SBATCH --time=01:00:00
#SBATCH --job-name=QC_analysis
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=qual_analysis_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Quality_control/fastqc/0.11.9

# run a FastQC analysis for each FASTQ file in argument 
for file in $@
do
	# run the QC
	# -o option : results are stored in the FASTQC subdirectory of the analysis directory
	fastqc -o $HOME/analysis/FASTQC $file
done

