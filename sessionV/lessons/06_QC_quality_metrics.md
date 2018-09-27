---
title: "ChIP Quality Assessment Metrics"
authors: "Mary Piper and Meeta Mistry"
date: "March 17, 2018"
---

Approximate time: 60 minutes

## Learning Objectives

* Discuss sources of low quality ChIP-seq data
* Understand strand cross-correlation
* Evaluate and interpret QC metrics and the cross-correlation plot

## ChIP-seq quality assessment

Prior to performing any downstream analyses with the results from a peak caller, it is best practice to assess the quality of your ChIP-seq data. What we are looking for is good quality ChIP enrichment over background.

<img src="../img/chip_workflow_june2017_step3.png" width="700">

## Strand cross-correlation

A very useful ChIP-seq quality metric that is independent of peak calling is strand cross-correlation. It is based on the fact that a high-quality ChIP-seq experiment will produce significant clustering of enriched DNA sequence tags at locations bound by the protein of interest, that present as a bimodal enrichment of reads on the forward and reverse strands.

The bimodal enrichment of reads is due to the following:

During the ChIP-seq experiment, the DNA is fragmented and the protein-bound fragments are immunoprecipitated. This generates DNA fragments containing the protein-bound region. 

The + strand of DNA is sequenced from the 5' end, generating the red reads in the figure below, and the - strand of DNA is sequenced from the 5' end, generating the blue reads in the figure below. 

<img src="../img/chip-fragments.png" width ="300">

*Nat Biotechnol. 2008 Dec; 26(12): 1351–1359*

Due to the sequencing of the 5' ends of the fragments, this results in an enrichment of reads from the + strand (blue in the image below) being slightly offset from the enrichment of reads from the - strand (red in the image below). We need to **determine the number of bases to shift the peaks to yield maximum correlation between the two peaks**, which **should** correspond to the predominant **fragment length**. We can calculate the shift yielding the maximum correlation using the **cross-correlation metric**.

**REPLACE THIS BIMODAL FIGURE**

### Computing Cross-Correlation Scores

The cross-correlation metric is computed as the **Pearson's linear correlation between the minus strand and the plus strand, after shifting minus strand by k base pairs.** Using a small genomic window as an example, let's walk through the details of the cross-correlation below.

**REPLACE ALL OF THESE FIGURES**

<img src="../img/cc1.png" width ="500">

<img src="../img/cc2.png" width ="500">

<img src="../img/cc3.png" width ="500">

In the end, we will have a table of values mapping each base pair shift to a Pearson correlation value. 

**Pearson correlation values are computed for *every peak* for *each chromosome* and values are multiplied by a scaling factor and then summed across all chromosomes.**

We can then **plot cross-correlation values (y-axis) against the shift value (x-axis)** to generate a cross-correlation plot.

The cross-correlation plot **typically produces two peaks**: a peak of enrichment corresponding to the predominant **fragment length** (highest correlation value) and a peak corresponding to the **read length** (“phantom” peak).

### Examples of Cross-Correlation Plots

**Strong signal**

High-quality ChIP-seq data sets tend to have a larger fragment-length peak compared with the read-length peak. An example of a **strong signal** is shown below using data from **CTCF (zinc-finger transcription factor)** in human cells. With a good antibody, transcription factors will typically result in 45,000 - 60,000 peaks. The red vertical line shows the dominant peak at the true peak shift, with a small bump at the blue vertical line representing the read length.

<img src="../img/ctcf.png" width="300"> 

**Weak signal**

An example of **weaker signal** is demonstrated below with a **Pol2** data. Here, this particular antibody is not very efficient and these are broad scattered peaks. We observe two peaks in the cross-correlation profile: one at the true peak shift (~185-200 bp) and the other at read length. For weak signal datasets, the **read-length peak will start to dominate**.

<img src="../img/Pol2.png" width ="300">

**No signal**

A failed experiment will resemble a cross-correlation plot using **input only**, in which we observe little or no peak for fragment length. Note in the example below the **strongest peak is the blue line (read length)** and there is basically no other significant peak in the profile. The absence of a peak is expected since there should be no significant clustering of fragments around specific target sites (except potentially weak biases in open chromatin regions depending on the protocol used).

<img src="../img/input.png" width="300"> 


### Quality Metrics Based on Cross-Correlation

Using the cross-correlation plot we can **compute metrics for assessing signal-to-noise ratios in a ChIP-seq experiment** and to ensure the fragment length is accurate based on the experimental design. Poor signal-to-noise and inaccurate fragment lengths can indicate problems with the ChIP-seq data. These metrics are described in more detail below:

#### Normalized strand cross-correlation coefficent (NSC):

The ratio of the maximal cross-correlation value divided by the background cross-correlation (minimum cross-correlation value over all possible strand shifts). 

<img src="https://latex.codecogs.com/gif.latex?\frac{max(CC&space;values)}{min(CCvalues)}" title="\frac{max(CC values)}{min(CCvalues)}" />

- higher NSC values indicate more enrichment (better signal:noise)
- low signal-to-noise: NSC values < 1.1
- minimum possible NSC value: 1 (no enrichment) 

#### Relative strand cross-correlation coefficient (RSC):

The ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. 

<img src="https://latex.codecogs.com/gif.latex?\frac{max(CCvalues)&space;-&space;background}{phantomCCvalue&space;-&space;background}" title="\frac{max(CCvalues) - background}{phantomCCvalue - background}" />

- high enrichment: RSC values > 1
- low signal-to-noise: RSC values < 0.8
- minimum possible RSC value: 0 (no enrichment)


> **NOTE:** Low NSC and RSC values can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth or a combination of these. Datasets with few binding sites (< 200) could be due to biological reasons, such as a factor that truly binds only a few sites in a particular tissue type, which would output low NSC and RSC scores.


## Additional Quality Metrics for ChIP-seq data

The [ENCODE consortium](https://genome.ucsc.edu/ENCODE/qualityMetrics.html) analyzes the quality of the data produced using a variety of metrics. We have already discussed metrics related to strand cross-correlation such as NSC and RSC. In this section, we will provide descriptions of additional metrics that **assess the distribution of signal within enriched regions, within/across expected annotations, across the whole genome, and within known artefact regions.**

> **NOTE**: For some of the metrics we give examples of what is considered a 'good measure' indicative of good quality data. Keep in mind that passing this threshold does not automatically mean that an experiment is successful and a values that fall below the threshold does not automatically mean failure!


### SSD

The SSD score is a measure used to indicate evidence of enrichment. It provides a measure of pileup across the genome and is computed by looking at the **standard deviation of signal pile-up along the genome normalised to the total number of reads**. An enriched sample typically has regions of significant pile-up so a higher SSD is more indicative of better enrichment. SSD scores are dependent on the degree of total genome wide signal pile-up and so are sensitive to regions of high signal found with Blacklisted regions as well as genuine ChIP enrichment. 


### FRiP: Fraction of reads in peaks

This value reports the percentage of reads that overlap within called peaks.  This is another good indication of how ”enriched” the sample is, or the success of the immunoprecipitation. It can be considered a **”signal-to-noise” measure of what proportion of the library consists of fragments from binding sites vs. background reads**. FRiP values will vary depending on the protein of interest. A typical good quality TF with successful enrichment would exhibit a FRiP around 5% or higher. A good quality PolII would exhibit a FRiP of 30% or higher. There are also known examples of	good data with FRiP < 1% (i.e. RNAPIII).
	

### Relative Enrichment of Genomic Intervals (REGI)

Using the genomic regions identified as called peaks, we can obtain **genomic annotation to show where reads map in terms of various genomic features**. We then evaluate the relative enrichment across these regions and make note of how this compares to what we expect for enrichment for our protein of interest.


### RiBL: Reads overlapping in Blacklisted Regions

It is important to keep track of and filter artifact regions that tend to show **artificially high signal** (excessive unstructured anomalous reads mapping). As such the DAC Blacklisted Regions track was generated for the ENCODE modENCODE consortia. The blacklisted regions **typically appear uniquely mappable so simple mappability filters do not remove them**. These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats. 

<img src="../img/blacklist.png" width="600">

These regions tend to have a very high ratio of multi-mapping to unique mapping reads and high variance in mappability. **The signal from blacklisted regions has been shown to contribute to confound peak callers and fragment length estimation.** The RiBL score then may act as a guide for the level of background signal in a ChIP or input and is found to be correlated with SSD in input samples and the read length cross coverage score in both input and ChIP samples. These regions represent around 0.5% of genome, yet can account for high proportion of total signal (> 10%).

> **NOTE:** If you had filtered out blacklisted regions before peak calling, and those filtered BAM files are used as input to `ChIPQC` you will not need to evaluate this metric.

> **How were the 'blacklists compiled?** These blacklists were empirically derived from large compendia of data using a combination of automated heuristics and manual curation. Blacklists were generated for various species including and genome versions including human, mouse, worm and fly. The lists can be [downloaded here.](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/). For human, they used 80 open chromatin tracks (DNase and FAIRE datasets) and 12 ChIP-seq input/control tracks spanning ~60 cell lines in total. These blacklists are applicable to functional genomic data based on short-read sequencing (20-100bp reads). These are not directly applicable to RNA-seq or any other transcriptome data types. 




***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
