#!/usr/bin/env bash

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# convert the gtf file of the meta-assembly into bed format (only 6 first columns)
# option -v : store the type of bed file (for TSS or polyAsite) given in argument
awk -v var=$1 'BEGIN{ OFS = "\t" } 
	
	# select only the transcripts in known chromosomal region (no patches/scaffolds...) for compatibility with TSS bed file and bedtools
	$1 ~ /chr/ {	

		# extract columns for the bed format and take a window of 100bp for the overlap
		if ($3 == "transcript") {
	
			# discriminate if plus or minus strand
			if ($7 == "+") {
				
				# take the 100bp window at the transcript start/end coordinate depending if overlap with TSS or polyAsite
				if (var == "TSS") {
					if ($4 > 50){
						print $1,$4-50,$4+50,substr($12, 2, length($12)-3),$6,$7
					}else{
						print $1,$4,$4+50,substr($12, 2, length($12)-3),$6,$7
					}
				} else if (var == "polyA") {
					print $1,$5-50,$5+50,substr($12, 2, length($12)-3),$6,$7
				}
			}else{
				# take the 100bp window from the opposite coordinates when considering the - strand
                                if (var == "TSS") {
                                        print $1,$5-50,$5+50,substr($12, 2, length($12)-3),$6,$7
                                } else if (var == "polyA") {
                                	if ($4 > 50){
                                                print $1,$4-50,$4+50,substr($12, 2, length($12)-3),$6,$7
                                        }else{
                                                print $1,$4,$4+50,substr($12, 2, length($12)-3),$6,$7
                                        }
				}
			}
		}
	}' $HOME/data/meta_assembly.gtf > $HOME/data/meta_assembly_$1.bed

