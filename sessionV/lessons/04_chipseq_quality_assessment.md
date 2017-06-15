---
title: "ChIP-Seq Quality Assessment"
author: "Mary Piper"
date: "June 12, 2017"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 1.5 hours

## Learning Objectives

* Generate enrichment and quality metrics for evaluating ChIP-Seq data
* Generate a report containing additional quality metrics and diagnostic plots

## ChIP-Seq quality assessment

Prior to performing any downstream analyses with the results from a peak caller, it is best practice to assess the quality of your ChIP-Seq data. What we are looking for is good	quality	ChIP-seq	enrichment	over	background. **To quantify the quality, we will be using two different tools: [ChIPQC](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html) and [phantompeakqualtools](https://code.google.com/archive/p/phantompeakqualtools/).**

### Sources of low quality ChiP-seq

* **The specifity of the antibody.** It is possible that there was poor reactivity against the intended target in the IP. Or perhaps the antibody is not specific enough, causing cross-reactivity with other DNA-associated proteins.
* **Degree of enrichment:** 
* **Biases during library preparation:** PCR amplification and fragment biases.

### Quality Metrics for ChIP-seq data

The [ENCODE consortium](https://genome.ucsc.edu/ENCODE/qualityMetrics.html) analyzes the quality of the data produced using a variety of metrics. In this section we will provide descriptions of what some of these metrics are, and what they are measuring. Then, we will introduce the tools to be able to compute these metrics on your own ChIP-seq data.

#### FRIP: Fraction of reads in peaks
A useful	first-pass to evaluate	the	success	of	the	immunoprecipita)on
• A good	quality	TF	>	5%	-- but also known	examples	of	good	data	with	FRiP	<	1%	RNAPIII	and	ZNF274

#### REGI: Relative enrichment in genomic intervals


#### Signal to Noise:

These metrics are useful in determining the strength of the signal relative to noise and to ensure the fragment length is accurate based on the experimental design. Poor signal-to-noise and inaccurate fragment lengths can indicate problems with the ChIP-Seq data. They are described in more detail below:

**Normalized strand cross-correlation coefficent (NSC): **

The ratio of the maximal cross-correlation value divided by the background cross-correlation (minimum cross-correlation value over all possible strand shifts). Higher values indicate more enrichment, values less than 1.1 are relatively low NSC scores, and the minimum possible value is 1 (no enrichment). Datasets with NSC values much less than 1.05 tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type or it could be due to poor quality).

**Relative strand cross-correlation coefficient (RSC):**

The ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. The minimum possible value is 0 (no signal), highly enriched experiments have values greater than 1, and values much less than 1 may indicate low quality. RSC values significantly low (< 0.8) tend to have low signal to noise and can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.






## `phantompeakqualtools` 

The [`phantompeakqualtools`](https://code.google.com/archive/p/phantompeakqualtools/) package is a tool used to compute enrichment and quality measures for ChIP-Seq data [[1](http://www.g3journal.org/content/4/2/209.full)]. We will be using the package to compute the predominant insert-size (fragment length) based on strand cross-correlation peak and data quality measures based on relative phantom peak.

### Set up

The `phantompeakqualtools` package is written as an R script, that uses `samtools` as a dependency. The package has various options that need to be specified when running from the command line. To get set up, we will need to start an interactive session, load the necessary modules and set up the directory structure:

```
$ bsub -Is -n 6 -q interactive bash

$ module load stats/R/3.2.1 seq/samtools/1.2

$ cd ~/ngs_course/chipseq/results

$ mkdir chip_qc

$ cd chip_qc
```

### Downloading `phantompeakqualtools`

To use this `phantompeakqualtools` package, we need to download it from the project website. On the [project website](https://code.google.com/archive/p/phantompeakqualtools/), click on the *Downloads* option on the left-hand side of the page. The *Downloads* page has all updates for the package, with the most recent being from 2013. 

Right-click on the link for the most recent update, and copy the link.

Download the `phantompeakqualtools` to your directory using `wget`:

```
$ wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz

$ ls
```

> **NOTE:** *You may be asked to choose a mirror. If so, just choose a location nearest to where you are located (i.e. in the northeast). Some mirrors can be slower than others for downloads depending on the server speed and distance to the server.*

You should see `ccQualityControl.v.1.1.tar.gz` appear in the folder. This is a compressed folder, to extract the contents we use the `tar -xzf` command:

```
$ tar -xzf ccQualityControl.v.1.1.tar.gz
```
The `tar` command offers a simple way to compress and uncompress entire directories. We are using the command to uncompress the `ccQualityControl.v.1.1.tar.gz` directory. 

The options included are:

`-x`: extract a tar archive (or tarball) file

`-z`: the file is a compressed gzip archive file

`-f`: file name of archive file (needs to precede the file name)

> **NOTE:** *To compress a directory, you would issue the same command, but replace -x with -c, which specifies to create a new tar archive (or tarball) file, and after the name of the tar file you would name the directory to be compressed*

You should now see a `phantompeakqualtools` folder. Let's explore the contents a bit:

```
$ cd phantompeakqualtools

$ ls -l
```
There should also be a `README.txt` which contains all the commands, options, and output descriptions. Let's check out the `README.txt`:

```
$ less README.txt
```
Note that there are two R scripts that are described in the README file. Both will compute the fragment length, and data quality characteristics based on cross-correlation analysis, but one is for use in situations where the duplicates have been removed (`run_spp_nodups.R`). This is the script we will be using.

### Installing R libraries

In the README you will have noticed an *INSTALLATION* section. We will need to install the R package, `caTools`, into our personal R library to run the script. To do this, first open up R:

```
$ R
```

Use the install.packages() function to install `caTools`:

```
> install.packages("caTools", lib="~/R/library")

# Choose a mirror near to your location (i.e. northeast). I chose the PA 1 mirror, which is number 117.

> quit()

```
 > **NOTE:** We do not need to install `spp` because the R module we have loaded has the package pre-installed.

### Running `phantompeakqualtools`

To obtain quality measures based on cross-correlation plots, we will be running the `run_spp_nodups.R` script from the command line which is a package built on SPP. This modified SPP package allows for determination of the cross-correlation peak and predominant fragment length in addition to peak calling. We will be using this package solely for obtaining these quality measures (no peak calling). 

The options that we will be using include:

* `-c`: full path and name (or URL) of tagAlign/BAM file
* `-savp`: save cross-correlation plot
* `-out`: will create and/or append to a file several important characteristics of the dataset described in more detail below.

```
## DO NOT RUN THIS
## THIS SCRIPT IS FOR COMUTING METRICS ON A SINGLE FILE
$ Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
```
>_**NOTE:** Even though the script is called `run_spp.R`, we aren't actually performing peak calling with SPP. 

From within the `phantompeakqualtools` directory, we will create output directories and use a 'for loop' to **run the script on every Nanog and Pouf51 BAM file**:

```
$ mkdir -p logs qual

$ for bam in ../../bowtie2/*Nanog*aln.bam ../../bowtie2/*Pou5f1*aln.bam
do 
bam2=`basename $bam _aln.bam`
Rscript run_spp_nodups.R -c=$bam -savp -out=qual/${bam2}.qual > logs/${bam2}.Rout
done
```

The for loop generates **three output files**. The **quality metrics** are written in a tab-delimited text file, and the **log files** contains the standard output text. A third file is created in the same directory as the BAM files. These are pdf files that contain the **cross-correlation** plot for each sample. Let's move those files into the appropriate output directory:

```
$ mv ../../bowtie2/*pdf qual  

```

To visualize the quality metrics (.qual) files more easily, we will concatenate the files together to create a single summary file that you can move over locally and open up with Excel.

```
$ cat qual/*qual > qual/phantompeaks_summary.xls
```
Let's use Filezilla or `scp` move the summary file over to our local machine for viewing.

### Quality metrics output

The qual files are tab-delimited with the columns containing the following information:

- COL1: Filename: tagAlign/BAM filename 
- COL2: numReads: effective sequencing depth i.e. total number of mapped reads in input file 
- COL3: estFragLen: comma separated strand cross-correlation peak(s) in decreasing order of correlation. (**NOTE:** The top 3 local maxima locations that are within 90% of the maximum cross-correlation value are output. In almost all cases, the top (first) value in the list represents the predominant fragment length.) 
- COL4: corr_estFragLen: comma separated strand cross-correlation value(s) in decreasing order (col2 follows the same order) 
- COL5: phantomPeak: Read length/phantom peak strand shift 
- COL6: corr_phantomPeak: Correlation value at phantom peak 
- COL7: argmin_corr: strand shift at which cross-correlation is lowest 
- COL8: min_corr: minimum value of cross-correlation 
- COL9: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8 
- COL10: Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8) 
- COL11: QualityTag: Quality tag based on thresholded RSC (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)

> **NOTE:** The most important metrics we are interested in are the values in columns 9 through 12, however these numbers are computed from values in the other columns.

All samples have quite high NSC values indicating more enrichment, a good signal to noise and a fair number of peaks. Nanog-rep2 has a comparably higher NSC value which might explain the increased number of peaks for that sample compared to the others. The RSC and quality tags further indicate good chip signal and a quality IP, yielding a very high quality tag. Based on these metrics, the samples look good for further analysis.

### Cross-correlation plots

The cross-correlation plots show the best estimate for strand shift and the cross-correlation values. This file can be viewed by transferring it to your local machine using FileZilla. Copy `H1hesc_Nanog_Rep1_chr12_aln.pdf` to your machine to view the strand shift. The cross correlation peak shows the highest cross-correlation at fragment length 105, **How does this compare to the one we generated using MACS?**.

<img src="../img/H1hesc_Nanog_Rep1_chr12_aln.png" width=400>


## `ChIPQC`: Additional quality metrics for ChIP-seq data




> **NOTE:** Many of the plots that were generated in the ChIPQC report can also be generated using [`deepTools`](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html), a suite of python tools developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq or MNase-seq. If you are interested in learning more we have a [lesson on quality assessment using deepTools](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionV/lessons/qc_deeptools.md).

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

