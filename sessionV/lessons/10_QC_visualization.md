---
title: "QC visualization of peaks with IGV"
author: "Meeta Mistry"
date: "Thursday July 29th, 2017"
---

Approximate time: 45 minutes

## Learning Objectives
* Generate bigWig files
* Use IGV to visualize BigWig, BED and data from ENCODE

## Visualization of ChIP-seq data

We are going to generate some bigWig files to visualize our data using IGV. Instead of using BAM files that can be large and cannot be normalized, we generate bigWig files that have been normalized for read count relative to the input control. We will be using the `bamCompare` command within [deepTools](https://deeptools.github.io/) for this: [instructions to generate normalized bigwigs with deepTools](https://github.com/fidelram/deepTools/wiki/Normalizations).

> "**The bigWig format** is useful for dense, continuous data that will be displayed in the Genome Browser as a graph. BigWig files are created from wiggle (wig) type files using the program wigToBigWig.
> 
> The bigWig files are in an indexed binary format. The main advantage of this format is that only those portions of the file needed to display a particular region are transferred to the Genome Browser server. Because of this, bigWig files have considerably faster display performance than regular wiggle files when working with large data sets. The bigWig file remains on your local web-accessible server (http, https or ftp), not on the UCSC server, and only the portion needed for the currently displayed chromosomal position is locally cached as a "sparse file"."
>
> -[https://genome.ucsc.edu/goldenpath/help/bigWig.html](https://genome.ucsc.edu/goldenpath/help/bigWig.html)

```bash
cd ~/ngs_course/chipseq/results/
mkdir visualization
```

```bash
module load seq/deeptools/2.4.0
bamCompare -h
```

```bash
bamCompare -b1 bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam -b2 bowtie2/H1hesc_Input_Rep1_chr12_aln.bam -o visualization/Nanog_Rep1_chr12.bw 2> visualization/Nanog_Rep1_bamcompare.log
bamCompare -b1 bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam -b2 bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -o visualization/Nanog_Rep2_chr12.bw 2> visualization/Nanog_Rep2_bamcompare.log
```

```bash
bamCompare -b1 bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam -b2 bowtie2/H1hesc_Input_Rep1_chr12_aln.bam -o visualization/Pou5f1_Rep1_chr12.bw 2> visualization/Pou5f1_Rep1_bamcompare.log
bamCompare -b1 bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam -b2 bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -o visualization/Pou5f1_Rep2_chr12.bw 2> visualization/Pou5f1_Rep2_bamcompare.log
```

Copy over the bigWig files to your laptop using filezilla or scp. Also copy over the BEDtools output files to your computer.

Start IGV and load the 2 rep1 files, and the overlap BED files. You will notice that there are positive and negative values on the track, what do you think this denotes in the context of normalization?

> You can generate a simple, non-normalized bigWig with `bamCoverage` and you won't see any negative values. 

Now load the `Nanog_vs_Pou5f1_edgeR_sig.bed` and `Nanog_vs_Pou5f1_deseq2_sig.bed` into IGV.

Finally, we are going to visually compare our output to the output from the full dataset from ENCODE, by loading that data from the IGV server.

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
