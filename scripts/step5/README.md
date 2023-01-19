# step5
* To create the sleuth objects necessary for the differential expression analysis with Sleuth (R package version 0.30.1, R version 4.2.1) based on the abundance results for each replicate from Kallisto, use the R script "create_sleuth_objects.R" in Rstudio or with R (compatibility problem on the cluster, cannot install R packages in an R session). \
It results in two sleuth objects, "so_transcripts.sleuth" and "so_genes.sleuth" in the data directory. 
* To do the DE analysis at the transcript level with Sleuth, use the R script "transcript_level_diff_analysis.R" in Rstudio. \
It results in two tables, "transcripts_DE.tsv" (transcript-level differential expression table) and "stat_transcripts_DE.tsv" (statistics about the DE), and one graph "volcano_plot_transcripts.jpeg" for visualisation in the analysis directory.
* To do the DE analysis at the gene-level with Sleuth, use the R script "gene_level_diff_analysis.R" in Rstudio. \
It results in two tables, "genes_DE.tsv" (transcript-level differential expression table) and "stat_genes_DE.tsv" (statistics about the DE), and two graphs, "Heat_map_genes.pdf" and "comparison_known_DE_genes.pdf" for visualisation in the analysis directory.

