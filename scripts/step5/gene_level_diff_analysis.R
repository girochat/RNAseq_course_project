library(sleuth)
library(tidyverse)
if (!require("cowplot", quietly = TRUE))
  install.packages("cowplot")
library(cowplot)

setwd("/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat")

# load the sleuth object built for the gene-level DE
so_G <- sleuth_load("data/so_genes.sleuth")


# Sleuth Test Results (Wald test)
#################################

# store the results of the Wald test and view the most significant DE genes
G_table_wt <- sleuth_results(so_G, 'conditionParaclonal', 'wt', 
                             show_all = FALSE)
G_significant_wt <- dplyr::filter(G_table_wt, qval < 0.05)
head(G_significant_wt)


# Statistics
#############

# create table of differential expression about the most significant genes
# pval : p-value from the Wald test
# qval : adjusted p-value with the Benjamini-Hochberg correction
# esti_log2FC : estimated log2 fold change 
DE_G <- data.frame(gene_name = G_significant_wt$target_id,
                   pval = G_significant_wt$pval,
                   qval = G_significant_wt$qval,
                   esti_log2FC = G_significant_wt$b)

# add the known status and single exon property to the table of DE
table_transcript_gene <- read.table("analysis/transcripts_genes_id.tsv.sorted", header = T)
DE_G <- merge(table_transcript_gene, DE_G)
DE_G <- unique(subset(DE_G, select = -transcript_id, MARGIN=1))

# save the table of differential expression in a tsv file
write_tsv(DE_G, "analysis/genes_DE.tsv", col_names = T)

# create table of statistics about the DE of transcripts
nb_DE <- c(length(G_significant_wt$target_id), sum(DE_G$known_status == 0))
nb_filtered <- c(length(G_table_wt$target_id), NA )

nb_up <- c(sum(DE_G$esti_log2FC > 1), 
           sum((DE_G$esti_log2FC > 1)&(DE_G$known_status == 0)))
nb_down <- c(sum(DE_G$esti_log2FC < -1), 
             sum((DE_G$esti_log2FC < -1)&(DE_G$known_status == 0)))
nb_single_exon <- c(sum(DE_G$single_exon == 1), 
                    sum((DE_G$single_exon == 1)&(DE_G$known_status == 0)))

stat_DE_G <- data.frame(Type = c("gene", "novel_gene"), 
                        filtered = nb_filtered,
                        sign_DE = nb_DE, single_exon = nb_single_exon,
                        up = nb_up, down = nb_down, 
                        row.names = 1)

# save the table of statistics
write_tsv(stat_DE_G, "analysis/stat_genes_DE.tsv")


# Visualisations 
################

# draw a heat map of the 30 most significant DE genes
most_significant_genes <- head(G_significant_wt, 30)
pdf("analysis/Heat_map_genes.pdf")
plot_transcript_heatmap(so_G, most_significant_genes$target_id, color_high='red',
                        color_low = 'blue', color_mid = 'white')
dev.off()


# function that draws the box plot of a specific gene DE across samples
draw_boxplot <- function(gene_name){
  return((plot_bootstrap(so_G, gene_name, units='scaled_reads_per_base', 
                 color_by = "condition") + labs(y = "Estimated counts")))
}

# display the DE of 12 known genes from previous experiments for comparison
plot <- apply(matrix(c("SOX2", "VIM", "MYCL", "EPCAM", "THY1", "CD44", "ZEB2", "SNAI2",
                       "CDH1", "IL6", "MMP2", "CD274"), nrow = 1), MARGIN = 2,
              FUN = draw_boxplot)

# store the 12 boxplots in a pdf
pdf("analysis/comparison_known_DE_genes.pdf")
title <- ggdraw() + draw_label("Comparisons - Gene DE", fontface='bold')
plot_grid(plot[1][[1]], plot[2][[1]], plot[3][[1]], plot[4][[1]], plot[5][[1]],
               plot[6][[1]], ncol = 2, nrow=3, labels = "AUTO")
plot_grid(plot[7][[1]], plot[8][[1]], plot[9][[1]], plot[10][[1]],
plot[11][[1]], plot[12][[1]], ncol = 2, nrow = 3, labels = "AUTO")
dev.off()


# Sleuth Test Results (LRT test)
################################## 

G_table_lrt <- sleuth_results(so_G, 'reduced:full', "lrt", show_all = F)
G_significant_lrt <- dplyr::filter(G_table_lrt, qval < 0.05)
head(G_significant_lrt)

bonferroni_adj_P <- G_table_lrt$pval*length(G_table_lrt$target_id)
G_table_lrt <- dplyr::mutate(G_table_lrt, adj_pval=bonferroni_adj_P)
G_significant_lrt <- dplyr::filter(G_table_lrt, adj_pval <= 0.05)
G_significant_lrt
