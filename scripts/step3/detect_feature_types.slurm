#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1M
#SBATCH --time=00:10:00
#SBATCH --job-name=detect_feature
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# detect what are the feature types in the gtf file
awk 'BEGIN{

	# initalise the array with a known feature type
	array_key["exon"]=1
        array_value[1]="exon"
        i=2

# select non-commented lines only
}!/^#/{
	# add unknown feature type in the array
        if (!($3 in array_key))
        {
                array_key[$3]=1
                array_value[i]=$3
                i++
        }
}END{
	# display all detected feature types
        for (element in array_value)
        {
                print array_value[element]
        }
}' $HOME/data/meta_assembly.gtf > $HOME/analysis/meta_assembly_features.txt

