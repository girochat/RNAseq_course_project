library(ggplot2)
library(gridExtra)

setwd("/Users/girochat/Documents/Giliane/Bioinformatics/Master/RNA\ Sequencing/Project")

# function to create dataframe with transcripts having a specific feature
create_df_feature <- function(table, col, feature, val, panel, exon = FALSE){
  if (exon == TRUE){
    df <- table[table$single_exon == 0 & table[, col] == val, c(1, col)]
  }else{
    df <- table[table[, col] == val, c(1, col)]
    
  }
  df[, 1] <- "Transcripts"
  df[, 2] <- feature
  df[, 3] <- panel
  colnames(df) <- c("Transcript", "feature", "panel")
  
  return(df)
}


# import the table of DE transcripts with all features
table <- read.delim("analysis/final_table_DE.tsv")

# Known transcripts
###################

# create dataframes of specific features to plot the bar chart with ggplot
# (only for the DE transcripts)
n <- length(table$transcript_id)
ALL <- data.frame(Transcript = rep("Transcripts", length=n),
                feature = rep("Total", length=n), 
                panel = rep("DE", length=n))

type_PC <- create_df_feature(table, 12, "Protein coding", "prot_coding", "DE")

type_lncRNA <- create_df_feature(table, 12, "lncRNA", "lncRNA", "DE")

type_other <- create_df_feature(table, 12, "Other", "other", "DE")

type_unknown <- create_df_feature(table, 12, "Unknown", "unknown", "DE")

# Novel transcripts
####################

# select only the DE transcripts that are novel
table_novel <- table[table$known_status == "Unknown", ]

# create dataframes of specific features to plot the bar chart with ggplot
# (only for the DE and novel transcripts)
novel_ALL <- create_df_feature(table, 3, "Total", "Unknown", "DE & Novel")

table_novel_TSS <- table_novel[table_novel$TSS == "Yes",]
novel_TSS_polyA <- create_df_feature(table_novel_TSS, 7, "TSS & polyA", "Yes", 
                                     panel = "DE & Novel")

novel_PCpot <- create_df_feature(table_novel, 8, "Potential prot. coding", "coding", 
                                     panel = "DE & Novel")

novel_lncRNApot <- create_df_feature(table_novel, 8, "Potential lncRNA", "noncoding", 
                                     panel = "DE & Novel")

novel_interG <- create_df_feature(table_novel, 10, "Intergenic", "Yes", 
                                  panel = "DE & Novel")


table_novel_TSS_polyA <- table_novel_TSS[table_novel_TSS$polyAsite == "Yes",]
novel_PC_TSS_polyA <- create_df_feature(table_novel_TSS_polyA, 8, "High conf. prot. coding", "coding", 
                                        panel = "DE & Novel")

novel_lncRNA_TSS_polyA <- create_df_feature(table_novel_TSS_polyA, 8, "High conf. lncRNA", "noncoding", 
                                        panel = "DE & Novel")


# merge all the resulting dataframe into one (both for known and unknown)
features_T <- rbind(ALL, type_PC, type_lncRNA, type_other, type_unknown)

features_novel_T <- rbind(novel_ALL, novel_TSS_polyA,novel_PCpot, novel_lncRNApot, 
                         novel_interG, novel_PC_TSS_polyA, novel_lncRNA_TSS_polyA)


# order the features as they should appear on the plot for the known transcripts
Features <- factor(features_T$feature, 
                   levels = c("Total", "Protein coding", "lncRNA", "Other",
                              "Unknown"))

# plot the bar chart for the known ones
plot_T <- ggplot(features_T) + 
  geom_bar(aes(x = Transcript, fill = Features), 
           position = position_dodge(preserve = 'single')) +
  theme_bw() + 
  labs(x = "") +
  ggtitle("DE")

# order the features as they should appear on the plot for the unknown transcripts
Feature <- factor(features_novel_T$feature, 
                   levels = c("Total", "TSS & polyA", "Intergenic",
                              "Potential prot. coding", "Potential lncRNA",
                              "High conf. prot. coding", "High conf. lncRNA"))

# plot the bar chart for the unknown transcripts
plot_novel_T <- ggplot(features_novel_T) + 
  geom_bar(aes(x = Transcript, fill = Feature),
           position = position_dodge(preserve = 'single')) +
  theme_bw() + 
  labs(x = "") +
  ggtitle("DE & Unknown")


# save the two-panels plot in jpeg format
jpeg("analysis/bar_chart.jpeg", height = 480, width = 900, quality = 90, 
     pointsize = 18)
grid.arrange(plot_T, plot_novel_T, ncol=2, nrow=1)
dev.off()

