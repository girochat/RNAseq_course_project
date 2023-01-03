#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=find_polyA
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/BEDTools/2.29.2

# find all transcripts that have an annotated polyA site in the downstream region
# option -u : outputs transcript that matches a polyA site only once
# option -s : compares transcript and polyA site on the same strand only  
bedtools intersect \
	-a $HOME/data/meta_assembly_polyA.bed \
	-b $HOME/data/polyAsite_human_pruned.bed \
	-u \
	-s > $HOME/analysis/transcripts_with_polyA.tsv
