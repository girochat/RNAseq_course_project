#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1M
#SBATCH --time=01:00:00
#SBATCH --job-name=count_exon
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# detect what are the feature types in the gtf file
awk 'BEGIN{
        array_key["exon"]=1
        array_value[1]="exon"
        i=2
}!/^#/{
        if (!($3 in array_key))
        {
                array_key[$3]=1
                array_value[i]=$3
                i++
        }
}END{
        for (element in array_value)
        {
                print array_value[element]
        }
}' $HOME/data/meta_assembly.gtf > $HOME/analysis/feature_types.txt

