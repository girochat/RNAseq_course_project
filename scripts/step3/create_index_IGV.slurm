#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=index_igv
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/samtools/1.10

samtools index -b $HOME/data/3_2.sorted.bam
