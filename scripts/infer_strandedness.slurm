#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --job-name=infer_strandedness
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Quality_control/RSeQC/4.0.0

# infer the strandedness of the library used in the project using RSeQC
# (RSeQC requires a reference genome in bed format)
infer_experiment.py -r $HOME/data/hg38_RefSeq.bed -i $HOME/data/P1.sorted.bam > $HOME/analysis/results_strandedness.txt
