#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=250M
#SBATCH --time=01:00:00
#SBATCH --job-name=QC_Analysis
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat/scripts/qual_analysis_%j.e

# add the software FastQC
module add UHTS/Quality_control/fastqc/0.11.9

# for each file given in argument to the script run a FastQC analysis and save the result in the analysis/ directory
for file in $@
do
	fastqc -o /data/courses/rnaseq_course/lncRNAs/Project2/users/grochat/analysis $file
done


