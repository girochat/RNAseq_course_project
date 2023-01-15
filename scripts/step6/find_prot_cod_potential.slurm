#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:01:00
#SBATCH --job-name=CPC2
#SBATCH --mail-user=giliane.rochat@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --error=error_%j.e

# path to the project directory
HOME=/data/courses/rnaseq_course/lncRNAs/Project2/users/grochat

# install CPC2 version 1.0.1 (python3 version)
cd $HOME
gzip -dc CPC2_standalone-1.0.1.tar.gz | tar xf -
cd CPC2_standalone-1.0.1
export CPC_HOME="$PWD"
cd libs/libsvm
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18
make clean && make

# run the protein coding potential analysis with CPC2
# option -r : analyse also the reverse strand 
cd $HOME/CPC2_standalone-1.0.1
bin/CPC2.py -r -i $HOME/data/meta_assembly.fa -o $HOME/analysis/transcripts_prot_coding.tsv

# remove CPC2 package
cd $HOME
rm -iR CPC2_standalone-1.0.1*

