#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1M
#SBATCH --time=00:15:00
#SBATCH --job-name=count_transcripts
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBACTH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# count the number of transcripts, genes and exons in the meta-assembly
cat $HOME/data/meta_assembly.gtf | awk 'BEGIN{
}!/^#/{

	# count the number of exons
        if ($3 == "exon")       
        {
                nb_exons++

		# count the number of transcripts mapped to a unique exon
		if ($14 == "\"1\";")
		{
			nb_unique_exons++
		} else if ($14 == "\"2\";")
		{
			nb_unique_exons--
		}
	
	# count the number of transcripts
        } else if ($3 == "transcript")
        {
                nb_transcripts++
		
		# count the number of known transcripts (with a GENCODE ID)
		gencode_prefix=substr($12, 2, 3)
		if (gencode_prefix == "ENS")
		{
			nb_known_transcripts++
		}
        }

	# count the number of genes
	if (!($10 in genes))
        {
                nb_genes++
                genes[$10]=1
        }
	
	# count the number of known genes (with a gene name)
        if (($13 == "gene_name")&&(!($14 in known_genes)))
        {	
                nb_known_genes++
        	known_genes[$14]=1
        }

}END{
	# count the number of novel genes and transcripts
	nb_novel_genes=(nb_genes-nb_known_genes)
	nb_novel_transcripts=(nb_transcripts-nb_known_transcripts)
       
	# display the statistics
	print "The meta-assembly contains :"
        print "\t",nb_genes,"genes"
	print "\t",nb_transcripts,"transcripts"
        print "\t",nb_exons,"exons","\n"
	print "In these, there are :"
        printf "\t" nb_novel_genes " novel genes (%.2f%%)\n",nb_novel_genes/nb_genes*100
        printf "\t" nb_novel_transcripts " novel transcripts (%.2f%%)\n",nb_novel_transcripts/nb_transcripts*100
	printf "\t" nb_unique_exons " transcripts with a unique exon (%.2f%%)\n",nb_unique_exons/nb_transcripts*100
}' > results_stat_assembly.txt


