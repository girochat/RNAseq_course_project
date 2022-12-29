#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name=detect_header
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# extract the file ID for the output file
file_ID=$(basename $1 .bed)

# get all possible headers (1st column) of the bed file given in argument
awk 'BEGIN{
	i=1
}{
	header=$1
	if (!(header in headers)){
		headers[header]=1
		headers_value[i]=header
		i++
	}
}END{
	print length(headers) 
	for (value in headers_value){
		print headers_value[value]
	} 
}' $1 > $HOME/data/${file_ID}_headers.txt
