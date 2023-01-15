#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=analyse_quant
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

for file in $@
do
	# extract the file ID for the output	
	file_ID=$(basename $file | awk -F "_final" '{print $1}')

	# determine the nb of genes, transcripts (novel and known) detected after quantification
	tail -n +2 $file | awk '{ 
		
		# do a quality check by summing up the tpm
		sum_tpm += $8

		nb_transcripts++
		nb_known_transcripts += $3
		if (!($2 in gene_name)) 
		{
			nb_genes++
			gene_name[$2] = 1
			if ($3 == 0)
			{
				nb_novel_genes += 1
			}
		}
		if (($3 == 0)&&($4 == 1))
		{
			nb_single_exon_T++
		}				
	}END{ 
		# get the number of novel transcripts
        	nb_novel_transcripts=(nb_transcripts-nb_known_transcripts)
       		
		# display the quality check of summing up the tpm
		print "The sum of the TPM values is : " sum_tpm,"\n"

        	# display the statistics
        	print "After the quantification of reads, there are :"
        	print "\t" nb_genes,"genes"
		printf "\t" nb_novel_genes " novel genes (%.2f%%)\n\n",nb_novel_genes/nb_genes*100
       		print "\t" nb_transcripts,"transcripts"
		printf "\t" nb_novel_transcripts " novel transcripts (%.2f%%)\n",nb_novel_transcripts/nb_transcripts*100
        	printf "\t" nb_single_exon_T " novel transcripts with a unique exon (%.2f%%)\n\n",nb_single_exon_T/nb_novel_transcripts*100
	}' > $HOME/analysis/abundance/${file_ID}/stat_abundance.txt
done


