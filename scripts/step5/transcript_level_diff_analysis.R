if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
if (!require("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
library('tidyverse')
library('gridExtra')

library(sleuth)

setwd("/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat")

# load the sleuth object built for the transcript-level DE
so <- sleuth_load("data/so_transcripts.sleuth")

# Sleuth Test Results (Wald test)
#################################

# store the results of the Wald test and view the most significant DE transcripts
T_table_wt <- sleuth_results(so, 'conditionParaclonal', 'wt', show_all = FALSE)
T_significant_wt <- dplyr::filter(T_table_wt, qval <= 0.05)
head(T_significant_wt) 

# Statistics
#############

# create table of differential expression about the most significant transcripts
# pval : p-value from the Wald test
# qval : adjusted p-value with the Benjamini-Hochberg correction
# esti_log2FC : estimated log2 fold change 
DE_T <- data.frame(transcript_id = T_significant_wt$target_id,
                     pval = T_significant_wt$pval,
                     qval = T_significant_wt$qval,
                     esti_log2FC = T_significant_wt$b)

# add the gene name, known status and single exon property to the table of DE
table_transcript_gene <- read.table("analysis/transcripts_genes_id.tsv.sorted", header = T)
DE_T <- merge(table_transcript_gene, DE_T)

# get the mean of the log2 of the observed counts (norm. and filtered) for 
# condition and control
matrix_T <- sleuth_to_matrix(so, "obs_norm", "est_counts")
matrix_T <- matrix_T[T_significant_wt$target_id, ]
mean_obs <- data.frame(transcript_id = T_significant_wt$target_id,
                       mean_obs_control = rowMeans(log2(matrix_T[, 4:6]+0.5)),
                       mean_obs_condition = rowMeans(log2(matrix_T[, 1:3]+0.5)))

# add the columns with the mean of the counts to the table
DE_T <- merge(DE_T, mean_obs)
DE_T <- DE_T[order(DE_T$pval), ]
head(DE_T)

# save the table of differential expression in a tsv file
write_tsv(DE_T, "analysis/transcripts_DE.tsv", col_names = T)

# create table of statistics about the DE transcripts
nb_DE <- c(length(T_significant_wt$target_id), sum(DE_T$known_status == 0))
nb_filtered <- c(length(T_table_wt$target_id), NA )

nb_up <- c(sum(DE_T$esti_log2FC > 1), 
           sum((DE_T$esti_log2FC > 1)&(DE_T$known_status == 0)))
nb_down <- c(sum(DE_T$esti_log2FC < -1), 
             sum((DE_T$esti_log2FC < -1)&(DE_T$known_status == 0)))
nb_single_exon <- c(sum(DE_T$single_exon == 1), 
                    sum((DE_T$single_exon == 1)&(DE_T$known_status == 0)))

stat_DE_T <- data.frame(Type = c("transcript", "novel_transcript"), 
                        filtered = nb_filtered,
                        sign_DE = nb_DE, single_exon = nb_single_exon,
                        up = nb_up, down = nb_down, 
                        row.names = 1)

# save the table of statistics
write_tsv(stat_DE_T, "analysis/stat_transcripts_DE.tsv")


# Visualisations
################

# get the gene name of the top 15 DE transcripts for the volcano plot
top <- data.frame(transcript_id=T_table_wt$target_id, qval = T_table_wt$qval)
top <- merge(top, table_transcript_gene)
top <- top[order(top$qval),]
top$gene_name[16:76977] <- ""
adj_qval <- -log10(T_table_wt$qval)

# draw a volcano plot of the DE transcripts
jpeg("analysis/volcano_plot_transcripts.jpeg")

volcano_plot <- plot_volcano(so, test = "conditionParaclonal", 
                             test_type='wt', sig_level = 0.05) + 
  theme_bw() +
  geom_text(aes(T_table_wt$b, adj_qval, label = top$gene_name), 
            hjust = 0, 
            vjust = 0,
            nudge_y = 0.8, 
            size=3, 
            inherit.aes = F,
            check_overlap = T) + 
  labs(x = "log2FC", title = "Volcano Plot", subtitle = "Transcript-level DE") +
  xlim(-15, 15) +
  geom_vline(aes(xintercept = 1), color = "black", linetype = "longdash") +
  geom_vline(aes(xintercept = -1), color = "black", linetype = "longdash") +
  geom_vline(aes(xintercept = 0), color = "lightgrey", linetype = "solid")

volcano_plot 

dev.off()


# Sleuth Test Results (LRT test)
#################################

# store the results of the LRT test and view the most significant DE transcripts
# with FDR cutoff of 5%
T_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
T_significant_lrt <- dplyr::filter(T_table_lrt, qval <= 0.05)
head(T_significant_lrt)


# do the bonferroni correction and visualise the most significant DE transcripts
bonferroni_adj_P <- T_table_lrt$pval*length(T_table_lrt$target_id)
T_table_lrt <- dplyr::mutate(T_table_lrt, adj_pval=bonferroni_adj_P)
T_significant_lrt <- dplyr::filter(T_table_lrt, adj_pval <= 0.05)
T_significant_lrt
