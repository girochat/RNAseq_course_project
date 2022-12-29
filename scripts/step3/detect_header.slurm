#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=detect_header
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

awk 'BEGIN{
	i=1
}{
	header=substr($1, 1, 3)
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
}' $HOME/data/meta_assembly.gtf > $HOME/analysis/assembly_headers.txt
