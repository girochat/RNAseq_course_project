#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# get the file ID for the output file
file_ID=$(basename $1 .bed)

# extract the 6 first columns of the bed file given in argument
awk 'BEGIN{ OFS = "\t" }{ print $1,$2,$3,$4,$5,$6 }' $1 > $HOME/data/${file_ID}_pruned.bed
