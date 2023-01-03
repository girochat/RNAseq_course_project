#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=find_ORF
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add SequenceAnalysis/GenePrediction/cpat/1.2.4

# find all protein-coding transcripts using CPAT
# option -x : hexamer table of the human genome provided by CPAT
# option -d : logit model of the human genome provided by CPAT
cpat.py -x $HOME/data/Human_Hexamer.tsv -d $HOME/data/Human_logitModel.RData -g $HOME/data/meta_assembly.bed -r $HOME/data/human_genome.fa -o $HOME/analysis/transcripts_ORF_bed

# use the cutoff recommended by CPAT to select only transcripts with significant protein coding potential
# cutoff : coding prob. >= 0.364
awk 'BEGIN{ print "transcript_id\tcoding_prob" }
{ 

	if ($4 >= 0.364)
	{
		print $1 "\t" $4
	}
}' $HOME/analysis/transcripts_ORF_bed.dat > $HOME/analysis/transcripts_ORF_bed.tsv

