#!/usr/bin/env bash

HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# select the columns with the transcript ID and the biotype from the reference annotation
awk !/^#/'{
	print substr($12, 2, length($12)-3) "\t" substr($20, 2, length($20)-3)
}' $HOME/data/reference_annotation.gtf > $HOME/analysis/transcripts_biotype.tsv

