#!/bin/bash

# USAGE: sh rnaseq_analysis_on_input_file.sh <fastq files> <number of cores> <path and name of output directory>
# This script will take a fastq file, number of cores (for STAR), and a path to the new analysis directory as input. It will perform the following steps on the fastq file and save all results to the specified directory.
	## starting with fastqc, 
	## followed by splice-aware alignment with STAR, 
	## generation of counts associated with genes using featureCounts.

# debugging with set -x [OPTIONAL]

# set -x

# assign the command line input to new variables
fq=$1
cores=$2
output_dir=$3

# shorten the name of the file

fname=$(basename $fq .fq)

echo "****Running rnaseq analysis on $fname****"

# Loading all the modules and adding featureCounts to the PATH

module load seq/STAR/2.5.3a
module load seq/fastqc/0.11.3

export PATH=/opt/bcbio/local/bin:$PATH

# make all of our output directories
	## The -p option means mkdir will create the whole path if it does not exist, and refrain from complaining if it does exist

mkdir -p $output_dir/trimmed_fastq $output_dir/STAR_alignment $output_dir/counts

# define variables to make modifications easier

# genome and gtf files that are likely to change
genome=/groups/hbctraining/ngs-data-analysis2016/rnaseq/reference_data/reference_STAR 
gtf=~/ngs_course/rnaseq/reference_data/chr1-hg19_genes.gtf

# output of alignment
align_out_prefix=$output_dir/STAR_alignment/${fname}_

# input and output of counts
bam_file=$output_dir/STAR_alignment/${fname}_Aligned.sortedByCoord.out.bam
counts=$output_dir/counts/${fname}.counts

# FastQC on trimmed file
echo "****FastQC on trimmed fastq****"
fastqc $fq

# Alignment with STAR

echo "****Running STAR alignment****"

STAR --runThreadN $cores \
--genomeDir $genome \
--readFilesIn $fq \
--outFileNamePrefix $align_out_prefix \
--outFilterMultimapNmax 10 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

# Counting reads with featureCounts

echo "****Running featureCounts****"

featureCounts -T $cores -s 2 -a $gtf -o $counts $bam_file
awk '{print $1"\t"$7}' $counts > $counts.txt





