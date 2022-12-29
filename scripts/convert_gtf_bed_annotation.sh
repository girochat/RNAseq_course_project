#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# convert the gtf file of the reference annotation into bed format (only 6 first columns)
awk 'BEGIN{ OFS = "\t" } 
	
	# select only the transcripts in known chromosomal region (no patches/scaffolds...) for compatibility with bedtools
	$1 ~ /chr/ {	
		# extract columns for the bed format 
		if ($3 == "transcript") {
			print $1,$4,$5,substr($12, 2, length($12)-3),$6,$7
		}
	}' $HOME/data/reference_annotation_ALL.gtf > $HOME/data/reference_annotation_ALL.bed
