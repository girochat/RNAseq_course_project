#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:15:00
#SBATCH --job-name=assembly
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Aligner/stringtie/1.3.3b

# perform a transcriptome assembly with StringTie for each BAM file in argument
for BAM_file in $@
do
	# extract ID of the BAM file for the output file
	output_name=$(basename ${BAM_file} | cut -d "." -f 1)
	
	# do the assembly with StringTie
	# option -o : output file
	# option --rf : first-strand stranded library
	# option -G : reference annotation to guide the assembly
	stringtie -o $HOME/data/${output_name}.gtf --rf -G $HOME/data/reference_annotation_ALL.gtf ${BAM_file}
done
