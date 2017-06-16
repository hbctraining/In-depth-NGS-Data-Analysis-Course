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


### `phantompeakqualtools`: quality metrics output

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


### `ChIPQC`: quality metrics report






> **NOTE:** Many of the plots that were generated in the ChIPQC report can also be generated using [`deepTools`](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html), a suite of python tools developed for the efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq or MNase-seq. If you are interested in learning more we have a [lesson on quality assessment using deepTools](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionV/lessons/qc_deeptools.md).

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

