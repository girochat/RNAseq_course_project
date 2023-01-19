# step2
* To create an index of the human reference genome GRCh38 (gencodegenes.org, release 21) with HISAT2 (version 2.2.1), use : \
`sbatch create_index_genome.slurm` \
The resulting index is stored in 8 different files ending with the suffix .ht2 in the references directory of the project 2 (`/data/courses/rnaseq/lncRNAs/Project2/references/`).
* To align the reads to the reference genome GRCh38 with HISAT2 (version 2.2.1) use the command: \
`sbatch launch_parallel_align.slurm ../../data/[3P]*R1*.fastq.gz ../../data/[3P]*R2*.fastq.gz`\
The result of the alignment is 6 SAM files, one for each replicate and a log file which can be found in the subdirectory "alignment" of the analysis directory. 
* To convert and sort the SAM file of each replicate into a BAM file with samtools (version 1.10), use : \
`sbatch convert_SAM_BAM.slurm ../../data/[3P]*.SAM` \
It results in 6 BAM files "3_2.sorted.bam", "3_4.sorted.bam", "3_7.sorted.bam", "P1.sorted.bam", "P2.sorted.bam", "P3.sorted.bam" in the data directory.
* To check the strandedness of the library used for sequencing with RSeQC (version 4.0.0) and the human reference genome hg38 RefSeq in bed format (sourceforge.net/projects/rseqc): \
`sbatch infer_strandedness.slurm` \
The result of the analysis can be found in the file "results_strandedness.txt" in the analysis directory.

