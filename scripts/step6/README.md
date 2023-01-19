# step6
* To verify the chromosome name format in any bed/gtf file (for compatibility using bedtools), use: \
`./detect_headers.sh ../../data/<file_id> bed/gtf` \
Any resulting file is stored in the same directory as the script with the suffix "_headers.txt".
* To prune any bed file (keeping only the first 6 columns), use: \
`./prune_bed_file.sh ../../data/<file_id>.bed` \
The corresponding pruned bed file (with the suffix "_pruned.bed") is stored in the data directory.
### TSS
* To convert the meta-assembly GTF file in BED format (6 columns) and select the appropriate coordinates for the overlap using bedtools depending on the search for TSS or polyAsite, use : \
`./convert_gtf_bed_assembly.sh TSS` \
`./convert_gtf_bed_assembly.sh polyA` \
The results are "meta_assembly_TSS.bed" and "meta_assembly_polyA.bed" in the data directory.
* To find transcripts that overlap known TSS regions using bedtools (version 2.29.2) and the TSS peaks file from FANTOM5, use: \
`sbatch find_TSS.slurm` \
The resulting table of transcripts "transcripts_TSS.tsv" can be found in the analysis directory.
### PolyAsite
* To modify the chromosome name format of the polyAsite bed file (format "chr#" needed for bedtools) and select only the 6 first columns, use:\
`./modify_headers_polyA_bed.sh`\
The modified bed file "polyAsite_human_pruned.bed" is stored in the data directory.
* To find transcripts that overlap known polyAsite regions using bedtools (version 2.29.2) and the polyAsite bed file from PolyASite, use : \
`sbatch find_polyA.slurm` \
The resulting table of transcripts "transcripts_polyA.tsv" can be found in the analysis directory.
### Protein coding potential
* To calculate the protein coding potential of the transcripts using CPC2 (version 1.0.1 - python3 version), use : \
`sbatch find_prot_cod_potential.slurm` \
It results in the table "transcripts_prot_coding.tsv" in the analysis directory.
### Intergenic regions
* To convert the gtf file of the reference annotation in bed format with 6 columns, use: \
`./convert_gtf_bed.sh ../../data/reference_annotation_ALL.gtf` \
The resulting bed file "reference_annotation_ALL.bed" can be found in the data directory.
* To find transcripts present only in intergenic genomic regions (not overlapping the reference annotation) using bedtools (version 2.29.2), use : \
`sbatch find_intergenic.slurm` \
The resulting table of transcripts "transcripts_intergenic.tsv" can be found in the analysis directory.
### lncRNA regions
* To convert the gtf file of the reference annotation for lncRNA regions in bed format with 6 columns, use : \
`./convert_gtf_bed.sh ../../data/reference_annotation_lncRNA.gtf` \
The resulting bed file "reference_annotation_lncRNA.bed" can be found in the data directory.
* To find transcripts present in regions that code for known lncRNAs using bedtools (version 2.29.2): \
`sbatch find_overlap_lncRNA.slurm` \
The resulting table of transcripts "transcripts_lncRNA_regions.tsv" can be found in the analysis directory.

### Integrative analysis
* To get the length of the transcripts along with their ID and gene name, use the script: \
`sbatch fetch_transcripts_length.slurm` \
The result "transcripts_length.tsv" is stored in the analysis directory.
* To get a table with the biotype of the transcripts from the reference annotation: \
`./get_biotype_ref_annot.sh` \
It results in a table of transcript ID and their corresponding biotype "transcripts_biotype.tsv" in the analysis directory.
* To create a final table summarising all possible features of a transcript, use the R script "create_final_table_transcripts.R" in Rstudio (version 4.2.1). To obtain statistics about the different analysed features, use the R script "analyse_features.R". \
The final tables for all transcripts "final_table_all.tsv" and for the DE ones "final_table_DE.tsv" are stored in the analysis directory, as well as the tables of statistics for the known ("final_stat_known_T.tsv") and unknown ("final_stat_unknown_T.tsv") transcripts.

