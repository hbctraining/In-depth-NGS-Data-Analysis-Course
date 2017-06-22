---
title: "ChIP-Seq Quality Assessment"
author: "Mary Piper, Meeta Mistry"
date: "June 12, 2017"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 1.5 hours

## Learning Objectives

* Discuss other quality metrics for evaluating ChIP-Seq data
* Generate a report containing quality metrics using `ChIPQC`
* Identify sources of low quality data



## Additional Quality Metrics for ChIP-seq data

The [ENCODE consortium](https://genome.ucsc.edu/ENCODE/qualityMetrics.html) analyzes the quality of the data produced using a variety of metrics. In this section we will provide descriptions of what some of these metrics are, and what they are measuring. Then, we will introduce the tools to be able to compute these metrics on your own ChIP-seq data.


### SSD

The SSD score, as implemented in htSeqTools, is another indication of evidence of enrichment. It is computed by
looking at the standard deviation of signal pile-up along the genome normalised to the total number of reads. An
enriched sample typically has regions of significant pile-up so a higher SSD is more indicative of better enrichment. SSD
scores are dependent on the degree of total genome wide signal pile-up and so are sensitive to regions of high signal found
with Blacklisted regions as well as genuine ChIP enrichment. Here the first CTCF replicate shows the highest SSD score
and so greatest enrichment for depth of signal.
The final two metrics report the percentage of reads in different regions of interest. The first 

### FRiP: Fraction of reads in peaks

This value reports the percentage of reads that overlap within called peaks.  This is another good indication of how ”enriched” the sample is, or the success of the immunoprecipitation. It can be considered a ”signal-to-noise” measure of what proportion of the library consists of fragments from binding sites vs. background reads. FRiP values for ChIPs around 5% or higher generally reflect a successful enrichment, but there are also known	examples	of	good	data	with	FRiP	<	1% (i.e. RNAPIII).

### RiBL: Reads overlappin blacklisted regions

The final measure reports the percentage of reads that overlapped blacklisted regions (RiBL). The signal from blacklisted
has been shown to contribute to confound peak callers and fragment length estimation as well as contribute to the read
length peak in cross coverage and cross coverage read length peak ([2]). The RiBL score then may act as a guide for
the level of background signal in a ChIP or input and is found to be correlated with SSD in input samples and the read
length cross coverage score in both input and ChIP samples







## `ChIPQC`: quality metrics report

### Running `ChIPQC`

> **NOTE:** This next section assumes you have the `ChIPQC` package (vChIPQC_1.10.3) installed for R 3.3.3. If you haven't done this please run the following lines of code before proceeding.
>
```
source("http://bioconductor.org/biocLite.R")
biocLite("ChIPQC")
```

1. Before you begin copy over the BAM files and the corresponding indices (`*.bam*`) from `/groups/hbctraining/ngs-data-analysis-longcourse/chipseq-trimmed/results/bowtie2` to your local laptop using `FileZilla`. 

2. Also, copy over your peak calls (`.narrowPeak`) from MACS2 for each file from `/groups/hbctraining/ngs-data-analysis-longcourse/chipseq-trimmed/results/macs2` to your local laptop.

*NOTE: students will be using alignment files and peak calls from in-class results in their HOME directories*.

Download the sample data sheet available from [this link](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/raw/may2017/sessionV/samplesheet_chr12.csv).

3. Open up RStudio. File --> 'New Project' --> New directory --> chipqc

4. Create directories for `data` and `meta`. In `data` create subdirectories for `bams` and `peakcalls`.

```
## Load libraries
library(ChIPQC)

## Load sample data
samples <- read.csv('meta/samplesheet_chr12.csv')

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19") 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: Nanog and Pou5f1", reportFolder="ChIPQCreport")

```
An example report can be found [here](https://u35207958.dl.dropboxusercontent.com/u/35207958/chipseq-devel/ChIPQCreport/ChIP%20QC%20report%3A%20Nanog%20and%20Pou5f1.html).

## Sources of low quality ChiP-seq

Once you have identified low quality samples, th next logical step is to troubleshoot what might be causing it.

**The specifity of the antibody.** 

The quality of a ChIP experiment is ultimately dictated by the specificity of the antibody and the degree of enrichment achieved in the affinity precipitation step [[1]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). Antibody deficiencies are of two main types:

1. Poor reactivity against the intended target
2. Non-specific antibody, causing cross-reactivity with other DNA-associated proteins

Antibodies directed against transcription factors must be characterized using both a primary (i.e immunoblot, immunofluorescence) and secondary characterization (i.e knockout of target protein, IP with multiple antibodies).

* **Biases during library preparation:** 

*PCR amplification:* Biases arise because DNA sequence content and length determine the kinetics of annealing and denaturing in each cycle of PCR. The combination of temperature profile, polymerase and buffer used during PCR can therefore lead to differential efficiencies in amplification between different sequences, which could be exacerbated with increasing PCR cycles. This is often manifest as a bias towards GC rich fragments [[2]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4473780/). **Limited use of PCR amplification is recommended because bias increases with every PCR cycle.**

*Fragment bias*: The way in which sonication is carried out can result in different fragment size distributions and, consequently, sample-specific chromatin configuration induced biases. As a result, it is not recommended to use a single input sample as a control for ChIP-seq peak calling if it is not sonicated together with the ChIP sample. 



***

> **NOTE:** Many of the plots that were generated in the ChIPQC report can also be generated using [`deepTools`](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html), a suite of python tools developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq or MNase-seq. If you are interested in learning more we have a [lesson on quality assessment using deepTools](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionV/lessons/qc_deeptools.md).

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

