#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name=convert_gtf
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBACTH --mail-type=end,fail
#SBATCH --error=error_%j.e

HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/samtools/1.10
module add UHTS/Assembler/cufflinks/2.2.1

# create an index of the reference genome for gffread
samtools faidx $HOME/data/human_genome.fa

# extract transcripts from transcriptome gtf file to fasta file
gffread -w $HOME/data/meta_assembly.fa -g $HOME/data/human_genome.fa $HOME/data/meta_assembly.gtf
 
