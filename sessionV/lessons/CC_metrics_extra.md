

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


ADD IMAGES ABOUT PHANTOM PEAK = READ LENGTH FROM KAYLEIGH

