#! /bin/bash    

# Usage: sh <name of script> <fastq file>
# This script accepts a trimmed fastq file on the command line and processes through the alignment workflow and returns a count file.

fq=$1

# location to genome reference FASTA file
     genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
     gtf=~/unix_oct2015/rnaseq_project/data/reference_data/chr1-hg19_genes.gtf

# make all of our output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist
     mkdir -p ~/unix_oct2015/rnaseq_project/results/STAR
     mkdir -p ~/unix_oct2015/rnaseq_project/results/counts

# set up our software environment...
    module load seq/STAR/2.4.0j
    module load seq/samtools/1.2
    module load seq/htseq/0.6.1p1

echo "Processing file $fq ..."

# grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)
    echo "basename is $base"


# set up output filenames and locations
    align_out=~/unix_oct2015/rnaseq_project/results/STAR/${base}_
    align_in=~/unix_oct2015/rnaseq_project/results/STAR/${base}_Aligned.sortedByCoord.out.bam
    counts=~/unix_oct2015/rnaseq_project/results/counts/${base}.counts

# Run STAR
STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

# Create BAM index
samtools index $align_in

# Count mapped reads
htseq-count --stranded reverse --format bam $align_in $gtf  >  $counts


