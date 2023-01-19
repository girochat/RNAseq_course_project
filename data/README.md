# Data Directory
This directory contains the raw data from the RNA sequencing process as well as any data used during or resulting from the steps of the project. 

#### Symbolic links to the fasta files
The files such as `<replicate>_R#.fastq.gz` are symbolic links to the provided raw data of the different replicates for the paraclonal and parental cell subpopulations on the IBU cluster.

#### Symbolic links to the human reference genome and its index files
The file `human_genome.fa` is a symbolic link to the human reference genome FASTA file "GRCh38 (regions = ALL)" obtained on gencodegenes.org, release 21 (10.2022).  
The files `index_genome.#.ht2` are symbolic links to the index files forming the index of the human genome "GRCh38" for the alignment with HISAT2.  
The file `human_genome.fa.fai` is the index of the human genome "GRCh38" made with samtools.

#### Reference files
* `hg38_RefSeq.bed.gz` was downloaded from sourceforge.net/projects/rseqc (10.2022), to check the strandedness of the library with RSeQC.
* `reference_annotation_ALL.gtf` is the reference annotation "Comprehensive gene annotation (regions = ALL)" available on genocodegenes.org, release 21 (11.2022).
* `reference_annotation_lncRNA.gtf` is the reference annotation "Long non-coding RNA gene annotation (regions = CHR)" available on genocodegenes.org, release 21 (01.2023).
* `TSS_hg38_liftover.bed` is the hg38 remapped BED file obtained from the "TSS_human.bed.gz" file available on the FANTOM5 website (https://fantom.gsc.riken.jp, (12.2022)) using the liftOver tool.
* `polyAsite_human.bed` : BED file with the polyA sites mapped on the GRCh38 human genome version, source : polyasite.unibas.ch (12.2022). 
* `Human_cutoff.txt`, `Human_Hexamer.tsv` and `Human_logitModel.RData` are the protein coding potential reference files provided by CPAT.

#### BAM files
These are the resulting files after aligning the reads of each replicate to the human genome.  
Format : `<replicate>.sorted.bam`

#### Assembly files
The GTF files `<replicate>.gtf` are the individual assemblies for each replicate. The merged meta-assembly is available in different formats (GTF, BED, FASTA).  
The file `meta_assembly.k_fai` is the index of the meta-assembly done with Kallisto for the quantification.

#### Sleuth objects
Two R objects were created with the Sleuth package in R for the differential analysis: one for the transcripts `so_transcripts.sleuth` and one for the genes `so_genes.sleuth`.

