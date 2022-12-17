#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1M
#SBATCH --time=00:15:00
#SBATCH --job-name=fetch_ID
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

cat $HOME/data/meta_assembly.gtf | awk '{
	
	# fetch each transcript ID and the corresponding gene name or ID
	if ($3 == "transcript")
	{	
		# extract the ID of the transcript
        	transcript_id=substr($12, 2, length($12)-3) "\t"
        	
		# extract the gene name if it exists otherwise the gene ID
		if ($13 == "gene_name"){
                	transcript_id = transcript_id substr($14, 2, length($14)-3) "\t" 1
        	}else{
        		transcript_id = transcript_id substr($10, 2, length($10)-3) "\t" 0
        	}
		print transcript_id
	}
}' > $HOME/analysis/transcripts_genes_id.tsv
