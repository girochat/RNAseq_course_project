setwd("/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat/")

# create a final table with all possible features found for the transcripts :
# DE results, TSS, polyAsite, intergenic, lncRNA regions, protein coding potential, known biotype 

# function to import table about a feature and add a column for the type of feature
import_table <- function(filename, col, feature_type){
  table <- read.delim(filename, header = F)[, col]
  
  # add column with 'Yes' to indicate to feature is present
  column <- rep("Yes", length(table))
  table <- data.frame(table, column)
  colnames(table) <- c("transcript_id", feature_type)
  return(table)
}


# import table with all transcripts found in the meta-assembly
table_T <- read.delim("analysis/transcripts_genes_id.tsv.sorted")
table_T[table_T$known_status == 0, 3] <- "Unknown"
table_T[table_T$known_status == 1, 3] <- "Known"
table_T[table_T$single_exon == 0, 4] <- "No"
table_T[table_T$single_exon == 1, 4] <- "Yes"

# import table about DE (keep only columns about DE)
table_DE <- (read.delim("analysis/transcripts_DE.tsv"))[, c(1, 5, 6, 7, 8, 9)]

# import table about TSS
table_TSS <- import_table("analysis/transcripts_TSS.tsv", 4, "TSS")

# import table about polyA
table_polyA <- import_table("analysis/transcripts_polyA.tsv", 4, "polyAsite")

# import table about protein coding potential
table_PC_pot <- (read.delim("analysis/transcripts_prot_cod.tsv"))[, c(1, 8)]
colnames(table_PC_pot) <- c("transcript_id", "PC_potential")

# import table about lncRNA potential (overlap with known lncRNA sequences)
table_lncRNA_pot <- import_table("analysis/transcripts_lncRNA_regions.tsv", 4, "lncRNA_pot")

# import table about intergenic
table_interG <- import_table("analysis/transcripts_intergenic.tsv", 4, "intergenic")

# import the table about the biotypes and select only the columns transcript ID and type
table_biotype <- unique((read.delim("analysis/transcripts_biotype.tsv"))[, c(1,2)])
colnames(table_biotype) <- c("transcript_id", "ref_biotype")

# merge all tables into a final one with NA or 'No' values when the feature is absent 
final_table <- merge(table_T, table_TSS, all.x = T)
final_table$TSS[is.na(final_table$TSS)] <- "No"

final_table <- merge(final_table, table_polyA, all.x = T)
final_table$polyAsite[is.na(final_table$polyAsite)] <- "No"

final_table <- merge(final_table, table_PC_pot, all.x = T)

final_table <- merge(final_table, table_lncRNA_pot, all.x = T)
final_table$lncRNA_pot[is.na(final_table$lncRNA_pot)] <- "No"

final_table <- merge(final_table, table_interG, all.x = T)
final_table$intergenic[is.na(final_table$intergenic)] <- "No"

final_table <- merge(final_table, table_biotype, all.x =T)
final_table[is.na(final_table$ref_biotype), 10] <- "unknown"
biotype <- unique(final_table$ref_biotype)

# add an additional column grouping the biotypes into four subgroups :
# prot_coding, lncRNA, other, unknown
final_table["type"] <- "other"
final_table[final_table$ref_biotype == "unknown", 11] <- "unknown"
final_table[final_table$ref_biotype == "protein_coding", 11] <- "prot_coding"
final_table[final_table$ref_biotype == "lincRNA" | 
              final_table$ref_biotype == "non-coding" |
              final_table$ref_biotype == "antisense" |
              final_table$ref_biotype == "retained_intron" |
              final_table$ref_biotype == "sense_intronic" |
              final_table$ref_biotype == "sense_overlapping" |
              final_table$ref_biotype == "3prime_overlapping_ncrna", 11] <- "lncRNA"

final_table <- merge(final_table, table_DE, all.x = T)

# save the final table for all transcripts and for DE ones  
write.table(final_table, "analysis/final_table_all.tsv", sep = "\t", quote = F,
            row.names = F)
write.table(final_table[!is.na(final_table$pval), ], "analysis/final_table_DE.tsv", 
            sep = "\t", quote = F, row.names = F )


