---
title: "Peak calling with MACS2"
author: "Meeta Mistry"
date: "June 12th, 2017"
---

Contributors: Meeta Mistry, Radhika Khetani

Approximate time: 90 minutes

## Learning Objectives

* Understand the different components of the MACS2 algorithm
* Learn to use MACS2 for peak calling
* Interpret results from MACS2

# Peak Calling

Peak calling, the next step in our workflow, is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment. 

<img src="../img/chip_workflow_june2017_step2.png" width=700>


What we observe from the alignment files is a strand asymmetry with read densities on the +/- strand, centered around the binding site. The 5' ends of the selected fragments will form groups on the positive- and negative-strand. The distributions of these groups are then assessed using statistical measures and compared against background (input or mock IP samples) to determine if the binding site is significant.

<img src="../img/chip-fragments.png" width="300" align="middle"></div>

There are various tools that are available for peak calling. One of the more commonly used peack callers is MACS2, and we will demonstrate it in this session. *Note that in this Session the term 'tag' and sequence 'read' are used interchangeably.*

> **NOTE:** [SPP](http://www.nature.com.ezp-prod1.hul.harvard.edu/nbt/journal/v26/n12/full/nbt.1508.html) is also very commonly used for narrow peak calling. While we will not be going through the steps for this peak caller in this Session, we do have [a lesson on SPP](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionV/lessons/peak_calling_spp.md) that we encourage you to browse through if you are interested in learning more.

## MACS2

A commonly used tool for identifying transcription factor binding sites is named [Model-based Analysis of ChIP-Seq (MACS)](https://github.com/taoliu/MACS). The [MACS algorithm](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. Although it was developed for the detection of transcription factor binding sites it is also suited for larger regions.

MACS improves the spatial resolution of binding sites through **combining the information of both sequencing tag position and orientation.** MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity. The MACS workflow is depicted below. In this lesson, we will describe the steps in more detail.

<img src="../img/macs_workflow.png" width=500>


### Removing redundancy

MACS provides different options for dealing with **duplicate tags** at the exact same location, that is tags with the **same coordination and the same strand**. The default is to keep a single read at each location. The `auto` option, which is very commonly used, tells MACS to calculate the maximum tags at the exact same location based on binomal distribution using 1e-5 as pvalue cutoff. Another alternative is to set the `all` option keeps every tags. If an `integer` is given, at most this number of tags will be kept at the same location. This redundancy is addressed for both the ChIP and input samples.


### Modeling the shift size

The tag density around a true binding site should show a **bimodal enrichment pattern**. MACS takes advantage of this bimodal pattern to empirically model the shifting size to better locate the precise binding sites.


To find paired peaks to **build the model**, MACS first scans the whole dataset searching for highly significant enriched regions. *This is done only using the ChIP sample!* Given a sonication size (`bandwidth`) and a high-confidence fold-enrichment (`mfold`), MACS slides two `bandwidth` windows across the genome to find regions with **tags more than `mfold` enriched relative to a random tag genome distribution**. 

<img src="../img/model_shift.png" width=500>

MACS randomly samples 1,000 of these high-quality peaks, separates their Watson and Crick tags, and aligns them by the midpoint between their Watson and Crick tag centers. The distance between the modes of the Watson and Crick peaks in the alignment is defined as 'd' and represents the estimated fragment length. MACS shifts all the tags by d/2 toward the 3' ends to the most likely protein-DNA interaction sites.


### Scaling libraries

For experiments in which sequence depth differs between input and treatment samples, MACS linearly scales the **total control tag count to be the same as the total ChIP tag count**. The default behaviour is for the larger sample to be scaled down.

### Effective genome length

To calculate calculate λBG from tag count, MAC2 requires the **effective genome size** or the size of the genome that is mappable. Mappability is related to the uniqueness of the k-mers at a  particular position the genome. Low-complexity and repetitive regions have low uniqueness, which means low mappability. Therefore we need to provide the effective genome length to correct for the loss of true signals in low-mappable regions.

<img src="../img/mappable.png" width=300>

The mappability or uniqueness influences the average mapped depth (i.e if the effective genome length is small, the proportion of reads that map will be small). As shown in the table below mappability improves with increased read length. When low-mappable regions (e.g. a ratio  <  0.25) are of interest, it might be better to include multiple mapped reads or use paired-end reads.

<img src="../img/map_table.png" width=500>

### Peak detection

For ChIP-Seq experiments, tag distribution along the genome can be modeled by a Poisson distribution. After MACS shifts every tag, it then slides 2d windows across the genome to find candidate peaks with a significant tag enrichment (default is p < 10e-5). This is a Poisson distribution p-value based on λ. The Poisson is a one parameter model, where the parameter **λ is the expected number of reads in that window**.

<img src="../img/peak_detection.png" width=300>

Instead of using a uniform λ estimated from the whole genome, MACS uses a dynamic parameter, λlocal, defined for each candidate peak. The lambda parameter is estimated from the control sample and is deduced by **taking the maximum value across various window sizes: λlocal = max(λBG, λ1k, λ5k, λ10k).** In this way lambda captures the influence of local biases, and is robust against occasional low tag counts at small local regions. Possible sources for these biases include local chromatin structure, DNA amplification and sequencing bias, and genome copy number variation.

<img src="../img/lambda.png" width=300>

Overlapping enriched peaks are merged, and each tag position is extended d bases from its center. The location with the highest fragment pileup, hereafter referred to as the summit, is predicted as the precise binding location. The ratio between the ChIP-Seq tag count and λlocal is reported as the fold enrichment.


### Estimation of false discovery rate

Each peak is considered an independent test and thus, when we encounter thousands of significant peaks detected in a sample we have a multiple testing problem. In MACSv1.4, the FDR was determined empirically by exchanging the ChIP and control samples. However, in MACS2, p-values are now corrected for multiple comparison using the **Benjamini-Hochberg correction**.


## Running MACS2

We will be using the newest version of this tool, MACS2. The underlying algorithm for peak calling remains the same as before, but it comes with some enhancements in functionality. 


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
* `-g`: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.

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

> **NOTE:** Relaxing the q-value does not behave as expected in this case since it is partially tied to peak widths. Ideally, if you relaxed the thresholds, you would simply get more peaks but with MACS2 relaxing thresholds also results in wider peaks.

Now that we have a feel for the different ways we can tweak our command, let's set up the command for our run on Nanog-rep1:

```
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
 	-f BAM -g 1.3e+8 \
	-n Nanog-rep1 \
	--outdir macs2
```

The tool is quite verbose so you should see lines of text being printed to the terminal, describing each step that is being carried out. If that runs successfully, go ahead and **run the same command on the remaining samples**:

	 $ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Nanog-rep2
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Pou5f1-rep1
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --outdir macs2 -n Pou5f1-rep2

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

We can also generate plots using the R script file that was output by MACS2. There is a `_model.R` script in the directory. Let's load the R module and run the R script in the command line using the `Rscript` command as demonstrated below:


	$ module load stats/R/3.2.1
	$ Rscript Nanog-rep1_model.r
	
Now you should see a pdf file in your current directory by the same name. Create the plots for each of the samples and move them over to your laptop using `Filezilla`. 

Open up the pdf file for Nanog-rep1. The first plot illustrates **the distance between the modes from which the shift size was determined**. 

<img src="../img/model-macs.png" width="400">

The second plot is the  **cross-correlation plot**. This is a graphical representation of the Pearson correlation of positive- and negative- strand tag densities, shifting the strands relative to each other by increasing distance. We will talk about this in more detail in the next lesson.



***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
