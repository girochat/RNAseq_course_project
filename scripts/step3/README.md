# step3
* To perform the transcriptome assembly of each replicate with StringTie (version 1.3.3b), use : \
`sbatch assemble_transcripts.slurm ../../data/[3P]*.bam` \
It results in the six GTF files "3_2.gtf", "3_4.gtf", "3_7.gtf", "P1.gtf", "P2.gtf", "P3.gtf", in the data directory. 
* To merge all replicate assemblies into a meta-assembly with StringTie (version 1.3.3b), use : \
`sbatch merge_assembly.slurm` \
The result is the "meta_assembly.gtf" file in the data directory.
* To detect all possible feature types in the meta-assembly, use: \
`sbatch detect_feature_types.slurm` \
The resulting file "meta_assembly_features.txt" is stored in the same directory as the script.
* To perform statistics about the number of exons, transcripts and genes in the meta-assembly, use : \
`sbatch count_transcripts_genes.slurm` \
The result is stored in the file "results_stat_assembly.txt" in the analysis directory. 
* To create a table composed of the transcript ID, gene name, known status and single-exon property of every transcript of the meta-assembly, use : \
`sbatch fetch_transcripts_id.slurm` \
The result is stored in the "transcripts_genes_id.tsv" file in the analysis directory. 


