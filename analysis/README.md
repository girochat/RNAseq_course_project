# Analysis Directory

This directory contains the different results obtained at each step of the workflow. Here are some explanations about the naming syntax of the files/directories:

* abundance subdirectory : \
This subdirectory contains the abundance tables for each replicate output by Kallisto at the quantification step.
* alignment subdirectory : \
It contains the resulting log files after aligning the paired-end reads to the human reference genome.
* FASTQC subdirectory :\
This subdirectory contains the html files and the zip folder resulting from the qualitative analysis done with FASTQC.
* IGV subdirectory : \
Here can be found the screenshots of the IGV browser when verifying results (BAM files, meta-assembly, TSS cage peaks...) in the genome viewer.

Any file of the kind "transcripts_\<feature\>.tsv" contains the list of transcripts with the specific feature (ex: with a known TSS or polyA site)  or the list of transcripts with additional information (ex: length or biotype of the transcripts).
