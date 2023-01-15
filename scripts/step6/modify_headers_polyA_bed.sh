#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# modify the headers (chromosome ID in the 1st column) of the polyA site bed file for compatibility and select the 6 first columns
awk 'BEGIN{ OFS = "\t" }
	$1 ~ /[1-9XY]/ {
		print "chr" $1,$2,$3,$4,$5,$6
	}' $HOME/data/polyAsite_human.bed > $HOME/data/polyAsite_human_pruned.bed

