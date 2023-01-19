setwd("/Users/girochat/Documents/Giliane/Bioinformatics/Master/RNA\ Sequencing/Project")

table <- read.delim("analysis/final_table_all.tsv")

# create a table of ranked unknown transcripts as candidate for being true lncRNAs
# conditions to be in the table : 
# overlap with TSS & polyA site ; 
# non coding potential (CPC2) ;
# differentially expressed
table_ranked_transcripts <- table[table$TSS == "Yes" & table$polyAsite == "Yes" &
                                  table$PC_potential == "noncoding" &
                                    table$known_status == "Unknown" &
                                    !is.na(table$pval) & 
                                    table$single_exon == "No", ]

# order the table from lowest to highest qval
table_ranked_transcripts <- table_ranked_transcripts[order(table_ranked_transcripts$qval), ]

# save the table
write.table(table_ranked_transcripts, "analysis/ranked_candidates.tsv", sep = "\t", quote = F,
            row.names = F)
