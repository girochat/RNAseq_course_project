# install RHDF5 package with BiocManager from bioconductor.org
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.16")
BiocManager::install("rhdf5")

# install devtools from bioconductor.org
BiocManager::install("devtools")
library(devtools)

# install sleuth package from pachterlab on github
devtools::install_github("pachterlab/sleuth")
library('sleuth')

# get sleuth version
packageVersion('sleuth')  # sleuth version 0.30.1

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
if (!require("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
library('tidyverse')
library('gridExtra')


set.seed("08122022")
setwd("/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat")

# first we store the metadata of our replicates in a dataframe for sleuth to identify 
# which files are to be compared (parental (control) with paraclonal)

# get the sample IDs from the data directory
sample_ID <- dir(file.path("analysis", "abundance"))
sample_ID 

# get the path to the kallisto abundance.h5 file for each sample
kal_dirs <- file.path("analysis", "abundance", sample_ID)
kal_dirs


# create dataframe with the metadata in three columns
s2c <- data.frame(sample=sample_ID, condition=rep(c("Paraclonal", "Control_Parental"), 
                                                  each=3))
# add the path to the kallisto results of each replicate to the column 'path'
s2c <- dplyr::mutate(s2c, path=kal_dirs)
s2c

# Transcript-level DE
#####################

# prepare the sleuth object by importing kallisto's results, normalising and filtering counts 
# option transform_fun_counts : log base 2 for the FC
so <- sleuth_prep(s2c, extra_bootstrap_summary=T, 
                  transform_fun_counts = function(x) log2(x + 0.5))

# generate model where the expression is assumed to be different depending on the condition
so <- sleuth_fit(so, ~condition, 'full')

# generate model where the expression is assumed to be the same across samples (not impacted
# by a condition)
so <- sleuth_fit(so, ~1, 'reduced')

# run the likelihood ratio test with the generated models to detect DE transcript
so <- sleuth_lrt(so, "reduced", "full")

# run the Wald test to get the DE between the two conditions, with reference to
# the control (parental)
so <- sleuth_wt(so, 'conditionParaclonal', 'full')

# save the sleuth object
sleuth_save(so, "data/so_transcripts.sleuth")

# Gene-level DE
################

# import the table containing the gene IDs mapped to the transcript IDs and their known status
setwd("/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat")
t2g <- read.table("analysis/transcripts_genes_id.tsv.sorted", header = T)
head(t2g)

# add the gene IDs to the metadata for the sleuth object
t2g <- dplyr::rename(t2g, target_id = transcript_id)

# prepare the sleuth object by importing kallisto's results, normalising and filtering counts 
# option transform_fun_counts : log base 2 for the FC
# option 'gene_mode = T' to have gene-based modeling, plotting...
# option aggregation_column = "gene_name" to aggregate transcripts gene-wise
so_t2g <- sleuth_prep(s2c, extra_bootstrap_summary=T, target_mapping = t2g, 
                      aggregation_column = "gene_name", gene_mode = T, 
                      transform_fun_counts = function(x) log2(x + 0.5))

# generate model where expression is assumed to be condition-dependent
so_t2g <- sleuth_fit(so_t2g, ~condition, 'full')

# generate model where expression is assumed to be independent from any condition
so_t2g <- sleuth_fit(so_t2g, ~1, 'reduced')

# run the likelihood ratio test with the generated models 
so_t2g <- sleuth_lrt(so_t2g, 'reduced', 'full')

# run the Wald test to get the DE between the two conditions (condition = paraclonal)
so_t2g <- sleuth_wt(so_t2g, 'conditionParaclonal', 'full')

# save the sleuth object
sleuth_save(so_t2g, "data/so_genes.sleuth")


