#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=kallisto_index
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

module add UHTS/Analysis/kallisto/0.46.0

# create an index for the transcriptome fasta file with kallisto
# (note : script must be launched from the same directory as the data)
kallisto index -i meta_assembly.k_fai meta_assembly.fa

