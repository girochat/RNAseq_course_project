#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --job-name=merge_assembly
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Aligner/stringtie/1.3.3b

# merge all the separate assemblies into a meta-assembly
stringtie --merge -G $HOME/data/reference_annotation_ALL.gtf -o $HOME/data/meta_assembly.gtf $HOME/data/[3P].gtf
