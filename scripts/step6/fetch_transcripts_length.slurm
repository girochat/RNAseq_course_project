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

# create a table listing all transcripts of the meta-assembly with their gene name/ID, known status and if it has a single exon
cat $HOME/data/meta_assembly.gtf | awk 'BEGIN{ header=1 }{
	
	# fetch each transcript ID and the corresponding gene name or ID
	if ($3 == "transcript")
	{	
		if (header == 1) {
			print "transcript_id\tgene_name\tknown_status\tlength\tsingle_exon"
			header=0
		} else {

			if (unique_exon == 1){
				print transcript_id "\t" 1
			} else {
				print transcript_id "\t" 0
			}
			unique_exon=0
		}
 
		# extract the ID of the transcript
        	transcript_id = substr($12, 2, length($12)-3) "\t"
        	
		# extract the gene name if it exists otherwise the gene ID; add the corresponding known status
		if ($13 == "gene_name"){
                	transcript_id = transcript_id substr($14, 2, length($14)-3) "\t" 1
        	}else{
        		transcript_id = transcript_id substr($10, 2, length($10)-3) "\t" 0
        	}

		# extract the length of the transcript
		transcript_id = transcript_id "\t" ($5-$4)	
	}
	
	# determine if the transcript has a single exon
        if ($3 == "exon")       
        {
                if ($14 == "\"1\";")
                {
                        unique_exon++
                } else if ($14 == "\"2\";")
                {
                        unique_exon--
                }
        }
}END{ if (unique_exon == 1){
                                print transcript_id "\t" 1
                        } else {
                                print transcript_id "\t" 0
                        }
}' > transcripts_genes_id.tsv
