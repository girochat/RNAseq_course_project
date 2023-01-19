# step4
* To convert the gtf file of the meta-assembly in FASTA format with the gffread utility of Cufflinks (version 2.2.1): \
`sbatch convert_gtf_fasta.slurm` \
The FASTA file "meta_assembly.fa" can be found in the data directory.
* To perform an index of the meta-assembly (FASTA file) with Kallisto (version 0.46.0), use the command: \
(note : the script must be launched from the same directory than the data used in the script) \
`sbatch create_index_kallisto.slurm` \
The resulting index is the "meta_assembly.k_fai" file in the data directory.
* To quantify the reads of each replicate with Kallisto (version 0.46.0) relatively to the meta-assembly, use : \
(note : the script must be launched from the same directory than the data used in the script) \
`sbatch quantify_reads.slurm [3P]*R1*.fastq.gz [3P]*R2*.fastq.gz` \
The result of the quantification process with Kallisto comprises 3 different files ("abundance.h5", "abundance.tsv" and "run_info.json") which can be found in the abundance subdirectory corresponding to each replicate (abundance/\<replicate_ID\>) in the analysis directory.
* To select only transcripts with non-zero abundance and add their corresponding gene name, known status and single-exon property, use : \
`./create_final_abundance_table.sh ../../analysis/abundance/[3P]*/abundance.tsv`  
The final table for the read abundance "<replicate_ID\>_final_abundance.tsv" can be found in the abundance subdirectory of the analysis directory for each replicate.
* To perform brief statistics about the number of detected transcripts and genes (novel and known) after quantification, use : \
`./analyse_abundance.sh ../../analysis/abundance/[3P]*/[3P]*_final_*.tsv` \
The resulting statistics are stored in the file "stat_abundance.txt" in the abundance subdirectory of the analysis directory for each replicate.

