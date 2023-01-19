# Transcriptome assembly - RNA sequencing Project

## Project description :

This project aims at assembling the transcriptome of two particular cell subpopulations of non-small cell lung cancer cells (paraclonal and parental cells). Paraclonal cells are one characterised subtype of cells that develop from lung cancer stem cells and which has been associated with treatment resistance. The purpose is to identify novel and known genes or lncRNAs which are significantly differentially expressed compared to the control (parental cells) in order to direct further scientific investigation into those genomic regions and potentially improve our knowledge of the mechanisms that lead to tumor formation and treatment resistance. For the project, the raw data from the RNA sequencing of the paraclonal and the parental cell subpopulations was made available on the IBU cluster in the `/data/courses/rnaseq/lncRNAs/fastq/` directory. There are three replicates for each cell subpopulation (3_2/3_4/3_7 for the paraclonal replicates and P1/P2/P3 for the parental replicates).

## Workflow of the project :

* Step 1 - Qualitative and Quantitative analysis
The first step of the workflow is to check the quality of the reads after the sequencing process and to quantify the total number of reads for each replicate.
* Step 2 - Alignment
The second step involves mapping the reads to the human reference genome.
* Step 3 - Transcriptome assembly
The third step consists of assembling the previously aligned reads into biologically relevant genomic sequences (transcripts, exons...). To do so, a reference annotation of the human genome is used. All assemblies obtained from each replicate, parental and paraclonal taken together, is then merged into one meta-assembly.
* Step 4 - Quantification
For the fourth step of the analysis, the abundance of reads per transcript relative to the meta-assembly obtained at step 3 will be estimated. 
* Step 5 - Differential expression analysis
The fifth step of the workflow consists in the differential analysis of gene and transcript expression between the two cell subpopulations, the parental cell subpopulation serving as control.
* Step 6 - Integrative analysis
The sixth and last step involves searching for transcription start site and polyA site in the vicinity of the detected transcripts from the assembly. This analysis concerns in particular novel transcripts and genes previously identified. Possible presence of those novel transcripts in intergenic regions or lncRNAs annotated regions as well as their protein coding potential is also determined.


## Data and Directories :

The project directory is organised in three subdirectories : data, scripts and analysis. \
The scripts directory contains all the scripts used at each step of the workflow. \
The data directory contains symbolic links to the raw data from the sequencing and any data used at the different steps of the workflow (reference annotation file, index files, fantom5 TSS peaks file...). Any data generated during the steps of the workflow was stored under the data directory.\
The analysis directory contains the results of any analysis performed on the data obtained at the different steps.

Path for the repository of this project on the IBU cluster : `/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat/`
Link for the repository of this project on GitHub.com : <https://github.com/girochat/RNAseq_course_project>

