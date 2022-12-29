#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# for each abundance file given in argument, get the transcripts with non-zero read abundance and add the corresponding gene name
for file in $@
do
	# extract the generic name of the file
	file_id=$(echo "$file" | awk -F "_ab" '{print $1}')
	file_id=$(basename $file_id)	

	# extract the transcripts of TPM > 0
	awk '$5 != "0" { print $0 }' $file > temp.tsv
	
	# sort the detected transcripts 
	# note : use same LOCALE for sort and join for compatibility
	LANG=en_EN sort -k 1,1 temp.tsv > temp.tsv.sorted

	# add the header for the detected transcripts
	sed '1i\
target_id	length	eff_length	est_counts	tpm' temp.tsv.sorted > temp.tsv	

	# add the corresponding gene name, its known status and its single-exon property to the detected transcript
	LANG=en_EN join --header -1 1 -2 1 $HOME/analysis/transcripts_genes_id.tsv.sorted temp.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8}' > $HOME/analysis/${file_id}_abundance/${file_id}_final_abundance.tsv

done 

rm temp.tsv temp.tsv.sorted
