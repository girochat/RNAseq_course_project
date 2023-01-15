#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --job-name=convert_SAM
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module load UHTS/Analysis/samtools/1.10 

# convert and sort each SAM file in argument with samtools
for file in $@
do
	# extract the input file ID for the output file
	name=$(basename $file | cut -d "." -f 1)
	
	# convert and sort the SAM file
	# option -o : specify output file
	samtools sort -o $HOME/data/${name}.sorted.bam $file
done
 
