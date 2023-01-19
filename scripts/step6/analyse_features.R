setwd("/Users/girochat/Documents/Giliane/Bioinformatics/Master/RNA\ Sequencing/Project")

# calculate statistics about the different analysed features
############################################################

# import the final table with all the analysed features
table <- read.delim("analysis/final_table_all.tsv")


# check how many known prot. coding and lncRNA were correctly identified by CPC2
amount_PC_T <- sum(table$type == "prot_coding")
amount_pot_PC_T <- sum(table$PC_potential == "coding" & table$type == "prot_coding")

amount_lncRNA_T <- sum(table$type == "lncRNA")
amount_pot_lncRNA_T <- sum(table$type == "lncRNA" & table$PC_potential == "noncoding")

round(amount_pot_PC_T / amount_PC_T * 100, 2)             
## => 69.92 % correctly identified protein coding transcripts

round(amount_pot_lncRNA_T / amount_lncRNA_T * 100, 2)     
## => 74.39 % correctly identified lncRNA


# Stat about known transcripts
##############################

# create column of statistics for the following features : TSS, polyA, TSS&polyA, protein coding, lncRNA 
# (stats are calculated for all the transcripts, the known ones, the known DE ones)
stat_TSS <- c(sum(table$TSS == "Yes"), 
                    sum(table$TSS == "Yes" & table$known_status == "Known"),
                    sum(table$TSS == "Yes" & !is.na(table$pval) & table$known_status == "Known"))

stat_polyA <- c(sum(table$polyAsite == "Yes"), 
                      sum(table$polyAsite == "Yes" & table$known_status == "Known"),
                      sum(table$polyAsite == "Yes" & !is.na(table$pval) & table$known_status == "Known"))

stat_TSS_polyA <- c(sum(table$TSS == "Yes" & table$polyAsite == "Yes"), 
                          sum(table$TSS == "Yes" & table$polyAsite == "Yes" & table$known_status == "Known"),
                          sum(table$TSS == "Yes" & table$polyAsite == "Yes" & !is.na(table$pval) & table$known_status == "Known"))

stat_prot_coding <- c(sum(table$type == "prot_coding"), 
                      sum(table$type == "prot_coding" & table$known_status == "Known"),
                      sum(table$type == "prot_coding" & table$known_status == "Known" & !is.na(table$pval)))

stat_lncRNA <- c(sum(table$type == "lncRNA"), 
                      sum(table$type == "lncRNA" & table$known_status == "Known"),
                      sum(table$type == "lncRNA" & table$known_status == "Known" & !is.na(table$pval)))

# determine the total amount of transcripts per row to calculate the %
nb_T_row <- c(length(table$transcript_id), 
                        sum(table$known_status == "Known"), 
                        sum(!is.na(table$pval) & table$known_status == "Known"))

# name the rows of the table of statistics
row_names <- c("Transcripts", 
               "Known Transcripts", 
               "Known DE transcripts", 
               "Transcripts (%)", 
               "Known transcripts (%)", 
               "Known DE transcripts (%)")

# create a dataframe with the columns of statistics
stat <- data.frame(Type = row_names,
                      TSS = c(stat_TSS, round(stat_TSS/nb_T_row*100, 2)),
                      polyA = c(stat_polyA, round(stat_polyA/nb_T_row*100, 2)),
                      TSS_polyA = c(stat_TSS_polyA, round(stat_TSS_polyA/nb_T_row*100, 2)),
                      prot_coding = c(stat_prot_coding, round(stat_prot_coding/nb_T_row*100, 2)),
                      lncRNA = c(stat_lncRNA, round(stat_lncRNA/nb_T_row*100, 2)))

# store the table of statistics
write.table(stat, "analysis/final_stat_known_T.tsv", sep = "\t", quote = F,
            row.names = F)


# Stat about unknown transcripts
################################

# select only the unknown transcripts
table_unknown <- table[table$known_status == "Unknown", ]

# verify if there are some DE transcripts of length < 200 (relevant for the lncRNAs search)
sum(table_unknown$length < 200 & !is.na(table_unknown$pval))   
# => result = 0

# create column of statistics for the following features : TSS&polyA, intergenic
# prot_cod_potential, lncRNA_potential, 
# TSS&polyA + prot_cod_pot, TSS&polyA + intergenic + prot_cod_pot,
# TSS&polyA + lncRNA_pot, TSS&polyA + intergenic + lncRNA_pot
# TSS&polyA + lncRNA_pot + lncRNAregions
# (stats are calculated for all the transcripts, the DE ones, the DE ones with a single exon)
stat_TSS_polyA_unknown <- c(sum(table_unknown$TSS == "Yes" & table_unknown$polyAsite == "Yes"),
                            sum(table_unknown$TSS == "Yes" & table_unknown$polyAsite == "Yes" & !is.na(table_unknown$pval)),
                            sum(table_unknown$TSS == "Yes" & table_unknown$polyAsite == "Yes" & !is.na(table_unknown$pval) & table_unknown$single_exon == "Yes"))

stat_intergenic_unknown <- c(sum(table_unknown$intergenic == "Yes"),
                             sum(table_unknown$intergenic == "Yes" & !is.na(table_unknown$pval)),
                             sum(table_unknown$intergenic == "Yes" & !is.na(table_unknown$pval) & table_unknown$single_exon == "Yes"))

stat_prot_cod_unknown <- c(sum(table_unknown$PC_potential == "coding"),
                           sum(table_unknown$PC_potential == "coding" & !is.na(table_unknown$pval)),
                           sum(table_unknown$PC_potential == "coding" & !is.na(table_unknown$pval) & table_unknown$single_exon == "Yes"))

stat_lncRNA_unknown <- c(sum(table_unknown$PC_potential == "noncoding"),
                         sum(table_unknown$PC_potential == "noncoding" & !is.na(table_unknown$pval)),
                         sum(table_unknown$PC_potential == "noncoding" & !is.na(table_unknown$pval) & table_unknown$single_exon == "Yes"))

# select the novel transcripts of high confidence (TSS+polyA) with a protein coding potential
table_high_conf_PC <- table_unknown[table_unknown$TSS == "Yes" & 
                                   table_unknown$polyAsite == "Yes" & 
                                   table_unknown$PC_potential == "coding", ]

stat_high_conf_PC <- c(length(table_high_conf_PC$transcript_id), 
                       sum(!is.na(table_high_conf_PC$pval)),
                       sum(!is.na(table_high_conf_PC$pval) & table_high_conf_PC$single_exon == "Yes"))

stat_high_conf_PC_intergenic <- c(sum(table_high_conf_PC$intergenic == "Yes"), 
                                  sum(table_high_conf_PC$intergenic == "Yes" & !is.na(table_high_conf_PC$pval)),
                                  sum(table_high_conf_PC$intergenic == "Yes" & !is.na(table_high_conf_PC$pval) & table_high_conf_PC$single_exon == "Yes"))

# select the novel transcripts of high confidence (TSS+polyA) with a lncRNA potential
table_high_conf_lncRNA <- table_unknown[table_unknown$TSS == "Yes" & 
                                        table_unknown$polyAsite == "Yes" & 
                                        table_unknown$PC_potential == "noncoding", ]

stat_high_conf_lncRNA <- c(length(table_high_conf_lncRNA$transcript_id), 
                           sum(!is.na(table_high_conf_lncRNA$pval)),
                           sum(!is.na(table_high_conf_lncRNA$pval) & table_high_conf_lncRNA$single_exon == "Yes"))

stat_high_conf_lncRNA_intergenic <- c(sum(table_high_conf_lncRNA$intergenic == "Yes"), 
                                      sum(table_high_conf_lncRNA$intergenic == "Yes" & !is.na(table_high_conf_lncRNA$pval)),
                                      sum(table_high_conf_lncRNA$intergenic == "Yes" & !is.na(table_high_conf_lncRNA$pval) & table_high_conf_lncRNA$single_exon == "Yes"))

stat_high_conf_lncRNA_lncRNAregions <- c(sum(table_high_conf_lncRNA$lncRNAregions == "Yes"), 
                                      sum(table_high_conf_lncRNA$lncRNAregions == "Yes" & !is.na(table_high_conf_lncRNA$pval)),
                                      sum(table_high_conf_lncRNA$lncRNAregions == "Yes" & !is.na(table_high_conf_lncRNA$pval) & table_high_conf_lncRNA$single_exon == "Yes"))


# determine the total amount of transcripts per row to calculate the %
nb_T_row_unknown <- c(length(table_unknown$transcript_id), 
                      sum(!is.na(table_unknown$pval)), 
                      sum(!is.na(table_unknown$pval) & table_unknown$single_exon == "Yes"))

# name the rows of the table of statistics
row_names_unknown <- c("Unknown transcripts", 
                       "DE unknown transcripts", 
                       "DE unknown transcripts with single exon", 
                       "Unknown transcripts (%)", 
                       "DE unknown transcripts (%)", 
                       "DE unknown transcripts with single exon (%)")

# create a dataframe with the columns of statistics
stat_unknown <- data.frame(Type = row_names_unknown,
                   TSS_polyA = c(stat_TSS_polyA_unknown, round(stat_TSS_polyA_unknown/nb_T_row_unknown*100, 2)),
                   intergenic = c(stat_intergenic_unknown, round(stat_intergenic_unknown/nb_T_row_unknown*100, 2)),
                   prot_cod_pot = c(stat_prot_cod_unknown, round(stat_prot_cod_unknown/nb_T_row_unknown*100, 2)),
                   lncRNA_pot = c(stat_lncRNA_unknown, round(stat_lncRNA_unknown/nb_T_row_unknown*100, 2)),
                   high_conf_PC = c(stat_high_conf_PC, round(stat_high_conf_PC/nb_T_row_unknown*100, 2)),
                   high_conf_PC_intergenic = c(stat_high_conf_PC_intergenic, round(stat_high_conf_PC_intergenic/nb_T_row_unknown*100, 2)), 
                   high_conf_lncRNA = c(stat_high_conf_lncRNA, round(stat_high_conf_lncRNA/nb_T_row_unknown*100, 2)),
                   high_conf_lncRNA_intergenic = c(stat_high_conf_lncRNA_intergenic, round(stat_high_conf_lncRNA_intergenic/nb_T_row_unknown*100, 2)),
                   high_conf_lncRNA_lncRNAregions = c(stat_high_conf_lncRNA_lncRNAregions, round(stat_high_conf_lncRNA_lncRNAregions/nb_T_row_unknown*100, 2)))

# store the table of statistics for the unknown transcripts
write.table(stat_unknown, "analysis/final_stat_unknown_T.tsv", sep = "\t", quote = F,
            row.names = F)
