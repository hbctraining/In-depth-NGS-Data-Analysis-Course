#! /bin/bash


# Make directory
mkdir -p results/STAR

# Run STAR
STAR --runThreadN 6 --genomeDir /groups/hbctraining/unix_oct2015_other/reference_STAR --readFilesIn data/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq  --outFileNamePrefix results/STAR/Mov10_oe_1_ --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

# Create BAM index
samtools index results/STAR/Mov10_oe_1__Aligned.sortedByCoord.out.bam

# Make directory
mkdir -p results/counts

# Count mapped reads
htseq-count --stranded reverse --format bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam data/reference_data/chr1-hg19_genes.gtf  >  results/counts/Mov10_oe_1.counts
