#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:30:00
#SBATCH --job-name=quantify_reads
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

module add UHTS/Analysis/kallisto/0.46.0

# store the read1/read2 pairs of files given in argument in an array
paired_files=($@)

# get the number of paired files
let nb_pairs=${#paired_files[@]}/2-1

for i in $(seq 0 $nb_pairs)
do
	# get the index to specify the second file of a pair
	let j=$i+${nb_pairs}+1
	
	# extract the generic name for the output files
	file_ID=$(basename ${paired_files[i]} _R1.fastq.gz)

	# quantify number of reads per transcript found in the assembly using Kallisto
	# option -i : index file
	# option --rf-stranded : specify library strandedness
	# option -o : output directory
	# option -b : number of bootstraps
	# (note : script must be launched from the same directory as the data)
	kallisto quant -i $HOME/data/meta_assembly.k_fai -o $HOME/analysis/${file_ID}_abundance -t 2 -b 100 --rf-stranded ${paired_files[i]} ${paired_files[j]}

done
