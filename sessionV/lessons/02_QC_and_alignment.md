---
title: "ChIP-Seq QC and alignment"
author: "Mary Piper, Radhika Khetani"
date: "June 11th, 2017"
---

Contributors: Mary Piper, Radhika Khetani

Approximate time: 1 hour

## Learning Objectives

* to use previous knowledge of quality control steps to perform FastQC and trimming
* to understand parameters and perform alignment using Bowtie2

# ChIP-Seq analysis 

Now that we have our files and directory structure, we are ready to begin our ChIP-Seq analysis. 

![workflow_QC](../img/chipseq_workflow_QC_partial.png)

## Quality Control
For any NGS analysis method, our first step is ensuring our reads are of good quality prior to aligning them to the reference genome. We will use FastQC to get a good idea of the overall quality of our data, to identify whether any samples appear to be outliers, to examine our data for contamination, and to determine a trimming strategy.

Let's run FastQC on all of our files. 

Start an interactive session with 4 cores if don't have one going, and change directories to the untrimmed_fastq folder.

```
$ cd ~/ngs_course/chipseq/raw_data 

$ module load seq/fastqc/0.11.3 

$ fastqc H1hesc_Input_Rep1_chr12.fastq 
```

Now, move all of the `fastqc` files to the `results/untrimmed_fastqc` directory:

`$ mv *fastqc* ../results/untrimmed_fastqc/`

Transfer the FastQC zip file for Input replicate 1 to your local machine using FileZilla and view the report.

![fastqc](../img/fastqc_input_rep1.png)

**ADD TRIMMOMATIC EXPLANATION HERE? SINCE WE NEVER COVERED THIS IN RNA-SEQ**

Based on the sequence quality plot, trimming should be performed from both ends of the sequences. We will use Trimmomatic to trim the reads from both ends of the sequence using the following parameters:

* `SE`: Single End reads
* `-threads`: number of threads / cores
* `-phread33`: quality score format
* `LEADING`: cut bases off the start of a read, if below a threshold quality
* `TRAILING`: cut bases off the end of a read, if below a threshold quality
* `MINLEN`: drop an entire read if it is below a specified length

Since we are only trimming a single file, we will run the command in the interactive session rather than creating a script:

```
$ module load seq/Trimmomatic/0.33

$ java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
-threads 4 \
-phred33 \
H1hesc_Input_Rep1_chr12.fastq \
../results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
LEADING:20 \
TRAILING:20 \
MINLEN:36
```

Let's see how much trimming improved our reads by running FastQC again:

`$ fastqc ../results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq`

Move the FastQC folders to the results directory for trimmed FastQC results:

`$ mv ../results/trimmed/*fastqc* ../results/trimmed_fastqc/`

Using Filezilla, transfer the file for the trimmed Input replicate 1 FastQC to your computer.

![trimmed_fastqc](../img/chipseq_trimmed_fastqc.png)

## Alignment

![workflow_align](../img/chipseq_workflow_align_partial.png)

Now that we have removed the poor quality sequences from our data, we are ready to align the reads to the reference genome. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) is a fast and accurate alignment tool that indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

*Bowtie2* supports gapped, local and paired-end alignment modes. It works best for reads that are at least 50 bp (shorter read lengths should use Bowtie1), and it can perform soft-clipping to remove poor quality bases [[1](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25)].

> _**NOTE:** Our reads are only 36 bp, so technically we should use Bowtie1. However, since it is rare that you will have sequencing reads with less than 50 bp, we will show you how to perform alignment using Bowtie2._

To perform peak calling for ChIP-Seq analysis, we need our alignment files to contain only **uniquely mapping reads** (no multi-mappers or duplicate reads) in order to increase confidence in site discovery and improve reproducibility. Since there is no parameter in Bowtie2 to keep only uniquely mapping reads, we will need to perform the following steps to generate alignment files containing only the uniquely mapping reads:

1. Create a bowtie2 index
2. Align reads with bowtie2 and output a SAM alignment file
3. Change alignment file format from SAM to BAM
4. Sort BAM file by read coordinate locations
5. Filter to keep only uniquely mapping reads (this will also remove any unmapped reads)

### Creating Bowtie2 index

To perform the Bowtie2 alignment, a genome index is required. **We previously generated the genome indexes for you**, and they exist in the `reference_data` directory.

However, if you needed to create a genome index yourself, you would use the following command:

```
# DO NOT RUN

bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>

# Can find indexes for the entire genome on Orchestra using following path: /groups/shared_databases/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/

```

### Aligning reads with Bowtie2

Since we have our indexes already created, we can get started with read alignment. Change directories to the `bowtie2` folder:

```

$ cd ~/ngs_course/chipseq/results/bowtie2

```

We will perform alignment on our single trimmed sample, `H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq`. Details on Bowtie2 and its functionality can be found in the [user manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using Bowtie2 are:

* `-p`: number of processors / cores
* `-q`: reads are in FASTQ format
* `-x`: /path/to/genome_indices_directory
* `-U`: /path/to/FASTQ_file
* `-S`: /path/to/output/SAM_file

```
$ bowtie2 -p 4 -q \
-x ~/ngs_course/chipseq/reference_data/chr12 \
-U ~/ngs_course/chipseq/results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
-S ~/ngs_course/chipseq/results/bowtie2/H1hesc_Input_Rep1_chr12_aln_unsorted.sam

```
> _**NOTE:** If you had added the bcbio path ino your `.bashrc` file last session you should be able to use bowtie2 without loading a module. If not, load the module using `module load seq/bowtie/2.2.4`_
>

### Changing file format from SAM to BAM

While the SAM alignment file output by Bowtie2 is human readable, we need a BAM alignment file for downstream tools. Therefore, we will use [Samtools](http://samtools.github.io) to convert the file formats. The command we will use is `samtools view` with the following parameters

* `-h`: include header in output
* `-S`: input is in SAM format
* `-b`: output BAM format
* `-o`: /path/to/output/file

```
$ samtools view -h -S \
-b H1hesc_Input_Rep1_chr12_aln_unsorted.sam \
-o H1hesc_Input_Rep1_chr12_aln_unsorted.bam
```

### Sorting BAM files by genomic coordinates

Before we can filter to keep the uniquely mapping reads, we need to sort our BAM alignment files by genomic coordinates. To perform this sort, we will use [Sambamba](http://lomereiter.github.io/sambamba/index.html), which is a tool that quickly processes BAM and SAM files. It is similar to SAMtools, but has unique functionality.
The command we will use is `sambamba sort` with the following parameters:

* `-t`: number of threads / cores
* `-o`: /path/to/output/file

```
$ sambamba sort -t 4 \
-o H1hesc_Input_Rep1_chr12_aln_sorted.bam \
H1hesc_Input_Rep1_chr12_aln_unsorted.bam 
```

### Filtering uniquely mapping reads

Finally, we can filter the uniquely mapped reads. We will use the `sambamba view` command with the following parameters:

* `-t`: number of threads / cores
* `-h`: print SAM header before reads
* `-f`: format of output file (default is SAM)
* `-F`: set [custom filter](https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax) - we will be using the filter to remove multimappers and unmapped reads.

```
$ sambamba view -h -t 4 -f bam \
-F "[XS] == null and not unmapped " H1hesc_Input_Rep1_chr12_aln_sorted.bam > H1hesc_Input_Rep1_chr12_aln.bam
```
We filtered out unmapped reads by specifying in the filter `not unmapped`. Also, among the reads that were aligned, we filtered out multimappers by specifying `[XS] == null`. 'XS' is a tag generated by Bowtie2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned and more than one alignment was found for the read [[1](http://computing.bio.cam.ac.uk/local/doc/bowtie2.html)].

Now that the alignment files contain only uniquely mapping reads, we are ready to perform peak calling.

> _**NOTE:** After performing read alignment, it's often useful to generate QC metrics for the alignment using tools such as [QualiMap](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#bam-qc) or [MultiQC](http://multiqc.info) prior to moving on to the next steps of the analysis._ 

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

