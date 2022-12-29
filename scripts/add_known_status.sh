#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# ignore the header
tail -n +2 $HOME/analysis/transcripts_genes_id.tsv.sorted > temp.tsv
	
# sort file by the transcript IDs 
# note : use same LOCALE for sort and join for compatibility
LANG=en_EN sort -k 4 $HOME/analysis/transcripts_with_TSS.tsv > $HOME/analysis/transcripts_with_TSS.tsv.sorted

# add the corresponding gene name and its known status to the transcripts with TSS
LANG=en_EN join -1 1 -2 4 temp.tsv $HOME/analysis/transcripts_with_TSS.tsv.sorted > $HOME/analysis/transcripts_TSS_known.tsv

rm temp.tsv

