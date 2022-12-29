#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

file_ID=$(basename $1 .gtf)

# convert the gtf file given in argument into bed format (only 6 first columns)
awk 'BEGIN{ OFS = "\t" } 
	
	# select only the transcripts in known chromosomal region (no patches/scaffolds...) for compatibility with bedtools
	$1 ~ /chr/ {	
	
		# extract columns for the bed format 
		if ($3 == "transcript") {
			print $1,$4,$5,substr($12, 2, length($12)-3),$6,$7
		}
	}' $1 > $HOME/data/${file_ID}.bed
