#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name=find_intergenic
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/BEDTools/2.29.2

# find transcripts that overlap intergenic genomic regions
# option -v : output transcripts that do NOT overlap 
# option -s : compare coordinates for sequences on the same strand only
bedtools intersect \
	-a $HOME/data/meta_assembly.bed \
	-b $HOME/data/reference_annotation_ALL.bed \
	-v \
	-s > $HOME/analysis/transcripts_intergenic.tsv

