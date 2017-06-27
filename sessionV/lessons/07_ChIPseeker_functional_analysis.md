---
title: "ChIPseeker for ChIP peak Annotation, Comparison, and Visualization"
author: "Meeta Mistry"
date: "June 12, 2017"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 1.5 hours

## Learning Objectives



## ChIPseeker

Now that we have a set of high confidence peaks for our samples, the next step is to **annotate our peaks to identify relative location relationship information between query peaks and genes/genomic features** to obtain some biological context. 

[ChIPseeker](http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) is an R package for annotating ChIP-seq data analysis. It supports annotating ChIP peaks and provides functions to visualize ChIP peaks coverage over chromosomes and profiles of peaks binding to TSS regions. Comparison of ChIP peak profiles and annotation are also supported, and can be useful to estimate how well biological replications are. Several visualization functions are implemented to visualize the peak annotation and statistical tools for enrichment analyses of functional annotations.


### Setting up 

1. Open up RStudio and open up the `chipseq-project` that we created previously.
2. Open up a new R script ('File' -> 'New File' -> 'Rscript'), and save it as `chipseeker.R`

> **NOTE:** This next section assumes you have the `ChIPseeker` package installed for R 3.3.3. You will also need one additional library for gene annotation. If you haven't done this please run the following lines of code before proceeding.
>
```
source("http://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
```

### Getting data 

As mentioned previously, these donwstream steps should be performed on your high confidence peak calls. While we have a set for our subsetted data, this set is rather small and will not result in anything meaningful in our functional analyses. **We have generated a set of high confidence peak calls using the full dataset.** These were obtained post-IDR analysis, (i.e. concordant peaks between replicates) and are provided in BED format which is optimal input for the ChIPseeker package. 

> **NOTE:** the number of peaks in these bed files are are significantly higher than what we observed with the subsetted data replicate analysis.

We will need to copy over the appropriate files from Orchestra to our laptop. You can do this using `FileZilla` or the `scp` command.

Move over the **BED files from Orchestra(`/groups/hbctraining/chip-seq/full-dataset/idr/*.bed`) to your laptop**. You will want to copy these files into your chipseq-project **into a new folder called `data/idr-bed`.**


Let's start by loading the libraries:

```r
# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
```

Now let's load all of the data. As input we need to provide the names of our BED files in a list format.

```r
# Load data
samplefiles <- list.files("data/idr-bed", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Nanog", "Pou5f1")

```

We need to assign annotation databases generated from UCSC to a variable:

	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

### Visualization

First, let's take a look at peak locations across the genome. The `covplot` function calculates **coverage of peak regions** across the genome and generates a figure to visualize this across chromosomes. We do this for the Nanog peaks and find a considerable number of peaks on all chromosomes.

```
# Assign peak data to variables
nanog <- readPeakFile(samplefiles[[1]])
pou5f1 <- readPeakFile(samplefiles[[2]])

# Plot covplot
covplot(nanog, weightCol="V5")

```

<img src="../img/covplot.png">


Using a window of +/- 1000bp around the TSS of genes we can plot the **density of read count frequency to see where binding is relative to the TSS** or each sample. This is similar to the plot in the ChIPQC report but with the flexibility to customize the plot a bit. We will plot both Nanog and Pou5f1 together to compare the two.

```
# Prepare the promotor regions
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

# Calculate the tag matrix
tagMatrixList <- lapply(as.list(samplefiles), getTagMatrix, windows=promoter)

## Profile plots
plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")
```
<img srce="../img/density_profileplots.png">

With these plots the confidence interval is estimated by bootstrap method (500 iterations) and is shown in the grey shading that follows each curve. The Nanog peaks exhibit a nice narrow peak at the TSS with small confidence intervals. Whereas the Pou5f1 peaks display a bit wider peak suggesting binding around the TSS with larger confidence intervals.

The **heatmap is another method of visualizing the read count frequency** relative to the TSS.

	# Plot heatmap
	tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)

<img src="../img/Rplot.png" width=500>

### Annotation


### Functional enrichment analysis







***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
