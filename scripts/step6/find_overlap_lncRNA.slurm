#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=find_lncRNAs
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/BEDTools/2.29.2

# find transcripts that overlap known lncRNA sequences
# option -s : compare coordinates for sequences on the same strand only
bedtools intersect \
        -a $HOME/data/meta_assembly.bed \
        -b $HOME/data/reference_annotation_lncRNA_2.bed \
        -u \
        -s > $HOME/analysis/transcripts_lncRNA_regions_2.tsv
