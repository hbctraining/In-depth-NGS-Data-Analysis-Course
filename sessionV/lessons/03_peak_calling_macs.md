---
title: "Peak calling with MACS2"
author: "Meeta Mistry"
date: "Tuesday, July 19th, 2016"
---

Contributors: Meeta Mistry, Radhika Khetani

Approximate time: 90 minutes

## Learning Objectives

* Understand the different components of the MACS2 algorithm
* Learn to use MACS2 for peak calling
* Interpret results from MACS2

## Peak Calling

Peak calling, the next step in our workflow, is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment. 

<div style="text-align:center"><img src="../img/workflow-peakcalling.png" width="400"></div>


From the alignment files (BAM), you typically observe reads/tags to be identified on each side of the binding site for the protein of interest. The 5' ends of the selected fragments will form groups on the positive- and negative-strand. The distributions of these groups are then assessed using statistical measures and compared against background (input or mock IP samples) to determine if the binding site is significant.


<div style="text-align:center"><img src="../img/chip-fragments.png" width="300" align="middle"></div>

There are various tools that are available for peak calling. Two of the popular ones we will demonstrate in this session are MACS2 and SPP. *Note that in this Session the term 'tag' and sequence 'read' are used interchangeably.*

## MACS2

A commonly used tool for identifying transcript factor binding sites is named [Model-based Analysis of ChIP-Seq (MACS)](https://github.com/taoliu/MACS). The [MACS algorithm](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity.


### Modeling the shift size
The tag density around a true binding site should show a bimodal enrichment pattern. MACS takes advantage of this bimodal pattern to empirically model the shifting size to better locate the precise binding sites.

![model](../img/model.png)

Given a sonication size (`bandwidth`) and a high-confidence fold-enrichment (`mfold`), MACS slides 2 bandwidth windows across the genome to find regions with tags more than mfold enriched relative to a random tag genome distribution. MACS randomly samples 1,000 of these high-quality peaks, separates their Watson and Crick tags, and aligns them by the midpoint between their Watson and Crick tag centers. The distance between the modes of the Watson and Crick peaks in the alignment is defined as 'd', and MACS shifts all the tags by d/2 toward the 3' ends to the most likely protein-DNA interaction sites.


### Peak detection

For experiments in which sequence depth differs between input and treatment samples, MACS linearly scales the total control tag count to be the same as the total ChIP tag count. Also, MACS allows each genomic position to contain no more than one tag and removes all the redundancies. 

To model the background noise, MACS uses a dynamic local Poisson distribution in which the lambda parameter is deduced by taking the maximum value λlocal = max(λBG, [λ1k,] λ5k, λ10k); this helps reduce the effects of local biases. Fold enrichment is then computed as the density of tags in a given peak compared to background λlocal parameter; this will help isolate the signal from the noise.


## Running MACS2

We will be using the newest version of this tool, MACS2. The underlying algorithm for peak calling remains the same as before, but it comes with some enhancements in functionality. One major difference on the computational end, is that for the older version (MACS) the FDR for peaks was empirically determined whereas with MACS2 the FDR is computed using the Benjamini-Hochberg method. ***Note: Relaxing the q-value does not behave as expected in this case since it is partially tied to peak widths. Ideally, if you relaxed the thresholds, you would simply get more peaks but with MACS2 relaxing thresholds also results in wider peaks.***


### Setting up

To run MACS2, we will first start an interactive session:

	$ bsub -Is -q interactive bash  
	
We will also need to create a directory for the output generated from MACS2:

	$ mkdir -p ~/ngs_course/chipseq/results/macs2
	
Now change directories to the `results` folder:

	$ cd ~/ngs_course/chipseq/results/
	
We only have the BAM file for our Input-rep1, but will need alignment information for **all 6 files**. We have generated the remaining BAM files for you, so **you will need to copy them over**:

	$ cp /groups/hbctraining/ngs-data-analysis-longcourse/chipseq/bowtie2/* bowtie2/

### MACS2 parameters

There are seven [major functions](https://github.com/taoliu/MACS#usage-of-macs2) available in MACS2 serving as sub-commands. We will only cover `callpeak` in this lesson, but if you can use `macs2 COMMAND -h` to find out more, if you are interested.

`callpeak` is the main function in MACS2 and can be invoked by typing `macs2 callpeak`. If you type this command without parameters, you will see a full description of commandline options. Here is a shorter list of the commonly used ones: 

**Input file options**

* `-t`: The IP data file (this is the only REQUIRED parameter for MACS)
* `-c`: The control or mock data file
* `-f`: format of input file; Default is "AUTO" which will allow MACS to decide the format automatically.
* `-g`: mappable genome size which is defined as the genome size which can be sequenced; some precompiled parameters provided.

**Output arguments**

* `--outdir`: MACS2 will save all output files into speficied folder for this option
* `-n`: The prefix string for output files
* `-B/--bdg`: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

**Shifting model arguments**

* `-s`: size of sequencing tags. Default, MACS will use the first 10 sequences from your input treatment file to determine it
* `--bw`: The bandwidth which is used to scan the genome ONLY for model building. Can be set to the expected sonication fragment size.
* `--mfold`: upper and lower limit for model building

**Peak calling arguments**

* `-q`: q-value (minimum FDR) cutoff
* `-p`: p-value cutoff (instead of q-value cutoff)
* `--nolambda`: do not consider the local bias/lambda at peak candidate regions
* `--broad`: broad peak calling


Now that we have a feel for the different ways we can tweak our command, let's set up the command for our run on Nanog-rep1:

```
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
 	-f BAM -g 1.3e+8 \
	--bdg --outdir macs2 \
	-n Nanog-rep1
```

The tool is quite verbose so you should see lines of text being printed to the terminal, describing each step that is being carried out. If that runs successfully, go ahead and **run the same command on the remaining samples**:

	 $ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Nanog-rep2
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Pou5f1-rep1
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Pou5f1-rep2

## MACS2 Output files

### File formats
Before we start exploring the output of MACS2, we'll briefly talk about some new file formats that we haven't yet encountered in this course.

**BED:**
The BED format consists of one line per feature, each containing 3-12 columns of data (whitespace-delimited or tab-delimited), plus optional track definition lines. This is a zero-based format.

The number of columns per line must be consistent throughout any single set of data in an annotation track. The first lines of of the file can consist of header lines. Header lines start with a hash (#) character, the word "browser," or the word "track."

The first three **required BED fields** are:

1. chrom - The name of the chromosome (e.g. chr3) or scaffold (e.g. scaffold10671)
2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
3. chromEnd - The ending position of the feature in the chromosome or scaffold. 

[Nine additional fields](http://useast.ensembl.org/info/website/upload/bed.html?redirect=no#optional) are optional.

**narrowPeak:**

A narrowPeak (.narrowPeak) file is used by the ENCODE project to provide called peaks of signal enrichement based on pooled, normalized (interpreted) data. It is a BED 6+4 format, which means the first 6 columns of a BED file with 4 additional fields:

7. signalValue - Measurement of overall enrichment for the region
8. pValue - Statistical significance (-log10)
9. qValue - Statistical significance using false discovery rate (-log10)
10. peak - Point-source called for this peak; 0-based offset from chromStart

**WIG format:**

Wiggle format (WIG) allows the display of continuous-valued data in a track format. Wiggle format is line-oriented. It is composed of declaration lines and data lines, and require a separate wiggle track definition line. There are two options for formatting wiggle data: variableStep and fixedStep. These formats were developed to allow the file to be written as compactly as possible.

**BedGraph format:**

The BedGraph format also allows display of continuous-valued data in track format. This display type is useful for probability scores and transcriptome data. This track type is similar to the wiggle (WIG) format, but unlike the wiggle format, data exported in the bedGraph format are preserved in their original state. For the purposes of visualization, these can be interchangeable.

### MACS2 output files

	$ cd macs2/
	
	$ ls -lh
	
There should be 6 files output to the results directory for each of the 4 samples, so a total of 24 files:

* `_peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
* `_peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
* `_summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended
* `_model.R`: an R script which you can use to produce a PDF image about the model based on your data and cross-correlation plot
* `_control_lambda.bdg`: bedGraph format for input sample
* `_treat_pileup.bdg`: bedGraph format for treatment sample

Let's first obtain a summary of how many peaks were called in each sample. We can do this by counting the lines in the `.narrowPeak` files:

	$ wc -l *.narrowPeak

We can also take a look at the plots, but first we will have to generate it. We can use the `Rscript` command to do this, but will first need to load the R module:

	$ module load stats/R/3.2.1
	$ Rscript Nanog-rep1_model.r
	
Now you should see a pdf file in your current directory by the same name. Create the plots for each of the samples and move them over to your laptop using `Filezilla`. 

Open up the pdf file for Nanog-rep1. The first plot illustrates the distance between the modes from which the shift size was determined. 

<img src="../img/model-macs.png" width="400">

The second plot is the cross-correlation plot, and we will be discussing that in a little more detail in the next (SPP) lesson.

**ADD CROSS-CORELATION DETAILS HERE IF WE TAKE OUT SPP**

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
