#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=7G
#SBATCH --time=00:30:00
#SBATCH --job-name=create_index
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --error=error_%j.e

$HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Aligner/hisat/2.2.1

# create the index of the reference genome with HISAT2
hisat2-build -p 8 $HOME/../../references/human_genome.fa $HOME/../../references/index_genome


