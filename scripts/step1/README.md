# step1  
* To perform a qualitative analysis of the FASTQ file of the parental and paraclonal cell replicates (both read1 and read2) with FASTQC (version 0.11.9), use : \
`sbatch check_quality_reads.slurm ../../data/[3P]*.fastq.gz` \
The results of the analysis can be viewed with the html file or in the zip archive file stored in the FASTQC subdirectory of the analysis directory for the pair of reads of each replicate.
* To count the number of reads found in each FASTQ file: \
`sbatch count_nb_reads.slurm ../../data/[3P]*.fastq.gz` \
The resulting table "results_nb_reads_fastq.tsv" can be found in the analysis directory.

