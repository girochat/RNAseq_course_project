#!/usr/bin/env bash

#SBATCH --array=0-5 %3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --job-name=align_reads
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# store paired files passed in argument in an array
# Note : read1 files (R1 in file path) are expected to come first, then read2 files (R2 in file path)
paired_files=($@)

# get the number of paired files to align
let nb_pairs=${#paired_files[@]}/2

# get the index to specify the second file of a pair
let i=${SLURM_ARRAY_TASK_ID}+${nb_pairs}

# launch the alignment for each pair of files
$HOME/scripts/align_to_ref_genome.sh ${paired_files[$SLURM_ARRAY_TASK_ID]} ${paired_files[$i]} 


