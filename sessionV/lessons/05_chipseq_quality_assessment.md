---
title: "ChIP-Seq alignment quality assessment"
author: "Mary Piper"
date: "Wednesday, July 20, 2016"
---

Contributors: Mary Piper

Approximate time: 1.5 hours

## Learning Objectives

* generate enrichment and quality measures for ChIP-Seq data
* assess the quality of alignments using coverage metrics and visualizations

# ChIP-Seq quality assessment

Prior to performing any analyses, it is best practice to assess the quality of your ChIP-Seq data for peak signal and alignment metrics. 

We will explore the *quality of the peaks* to determine the strength of the signal relative to noise and to ensure the fragment length is accurate based on the experimental design. Poor signal-to-noise and inaccurate fragment lengths can indicate problems with the ChIP-Seq data. 

For the *alignment quality*, we will investigate the read coverages for each sample and determine the variability in coverage per sample group. Replicate samples that vary greatly in where the reads stack up is indicative of a weak ChIP-Seq experiment. In addition, we can identify outlier samples or batch effects.

## Obtaining quality metrics using *phantompeakqualtools*

The *[phantompeakqualtools](https://code.google.com/archive/p/phantompeakqualtools/)* package allows for the generation of enrichment and quality measures for ChIP-Seq data [[1](http://www.g3journal.org/content/4/2/209.full)]. We will be using the package to compute the predominant insert-size (fragment length) based on strand cross-correlation peak and data quality measures based on relative phantom peak.

### Set up

The *phantompeakqualtools* package is written as an R script, `run_spp.R` that uses `samtools` as a dependency. The package has various options that can be specified when running from the command line. To get set up, we will need to start an interactive session, load the necessary modules and set up the directory structure:

```
$ bsub -Is -n 6 -q interactive bash

$ module load stats/R/3.2.1 seq/samtools/1.2

$ cd ~/ngs_course/chipseq/results

$ mkdir chip_qc

$ cd chip_qc
```
### Downloading *phantompeakqualtools*

To use this *phantompeakqualtools* package, we need to download it from the project website. On the [project website](https://code.google.com/archive/p/phantompeakqualtools/), click on the *Downloads* option on the left-hand side of the page. The *Downloads* page has all updates for the package, with the most recent being from 2013. 

Right-click on the link for the most recent update, and copy the link.

Download the *phantompeakqualtools* to your directory using `wget`:

```
$ wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz

$ ls
```

> ***NOTE:*** *You may be asked to choose a mirror. If so, just choose a location nearest to where you are located (i.e. in the northeast). Some mirrors can be slower than others for downloads depending on the server speed and distance to the server.*

You should see `ccQualityControl.v.1.1.tar.gz` appear in the folder. This is a compressed folder, to extract the contents we use the `tar -xzf` command:

```
$ tar -xzf ccQualityControl.v.1.1.tar.gz
```
The `tar` command offers a simple way to compress and uncompress entire directories. We are using the command to uncompress the `ccQualityControl.v.1.1.tar.gz` directory. 

The options included are:

`-x`: extract a tar archive (or tarball) file

`-z`: the file is a compressed gzip archive file

`-f`: file name of archive file (needs to precede the file name)

> ***NOTE:*** *To compress a directory, you would issue the same command, but replace -x with -c, which specifies to create a new tar archive (or tarball) file, and after the name of the tar file you would name the directory to be compressed*

You should now see a `phantompeakqualtools` folder. Let's explore the contents a bit:

```
$ cd phantompeakqualtools

$ ls -l
```

Note the script for generating the quality metrics, `run_spp.R`. There should also be a `README.txt` which contains all the commands, options, and output descriptions. Let's check out the `README.txt`:

```
$ less README.txt
```

### Installing R libraries

We will need to install the R package, `caTools`, into our personal R library to run the script:

```
$ R
```

In R, use the install.packages() function to install `caTools`:

```
> install.packages("caTools", lib="~/R/library")

# Choose a mirror near to your location (i.e. northeast). I chose the PA 1 mirror, which is number 117.

> quit()

```


### Running *phantompeakqualtools*

To obtain quality measures based on cross-correlation plots, we will be running the `run_spp.R` script from the command line which is a package built on SPP. This modified SPP package allows for determination of the cross-correlation peak and predominant fragment length without having to perform peak calling. We will be using this package solely for obtaining these quality measures. 

The options that we will be using include:

* `-c`: full path and name (or URL) of tagAlign/BAM file
* `-savp`: save cross-correlation plot
* `-out`: will create and/or append to a file several important characteristics of the dataset described in more detail below.

```
## DO NOT RUN THIS
$ Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
```
>_**NOTE:** Even though the script is called `run_spp.R`, we aren't actually performing peak calling with SPP. In addition, we could have run `run_spp_nodups.R` since we have already removed duplicates from our data and it would run slightly faster._

From the `phantompeakqualtools` directory, create output directories and run a 'for loop' to run the script on every Nanog and Pouf51 BAM file:

```
$ mkdir -p logs qual

$ for bam in ../../bowtie2/*Nanog*aln.bam ../../bowtie2/*Pou5f1*aln.bam
do 
bam2=`basename $bam _aln.bam`
Rscript run_spp.R -c=$bam -savp -out=qual/${bam2}.qual > logs/${bam2}.Rout
done
```

Now that we have our files created, we will move up a directory and organize our logs and quality output folders:

```
$ mv ../../bowtie2/*pdf qual  

# the pdf file gets created in the same directory as the input bam file, so we need to move it over.
```

To visualize the quality results (.qual) files more easily, we will concatenate the files together to create a single summary file that you can move over locally and open up with Excel.

```
$ cat qual/*qual > qual/phantompeaks_summary.qual
```
Let's use Filezilla or `scp` move the summary file over to our local machine for viewing.

#### Description of the quality information

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

Three of the more important values to observe are the NSC, RSC and QualityTag values:

**NSC:** values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type OR it could be due to poor quality)

**RSC:** values range from 0 to larger positive values. 1 is the critical threshold. RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.

**QualityTag:** A quick check of RSC, with negative values indicating poor signal to noise.

#### Cross-correlation plots

The cross-correlation plots show the best estimate for strand shift and the cross-correlation values. This file can be viewed by transferring it to your local machine using FileZilla. Copy `H1hesc_Nanog_Rep1_chr12_aln.pdf` to your machine to view the strand shift.

![phantom_peak](../img/H1hesc_Nanog_Rep1_chr12_aln.png)

The cross correlation peak shows the highest cross-correlation at fragment length 105, similar to what we found running spp. The Nanog rep1 has a high NSC value (greater than 1.1), indicating that Nanog rep1 exhibits good signal to noise and a fair number of peaks. The RSC and quality tags further indicate good chip signal and a quality IP, yielding a very high quality tag (>2). Based on these metrics, Nanog rep1 looks good for further analysis.

## Quality assessment using *deepTools*

Using the *[deepTools](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html)* suite of tools, we can assess the quality of our alignments for each of our samples by using several metrics to explore sample signal strength and read coverage similarity within and between experimental conditions.

Assessing and visualizing alignment quality using *deepTools* requires three steps: 

1. Indexing the BAM alignment files
2. Calculation of the read coverage scores using the `multiBamSummary` tool
3. Visualizing how read coverage scores compare between samples

To use deepTools, we will need an index (`.bai` file) for each of our BAM files using the `samtools index` tool. This has already been done for you, if you look inside your `bowtie2` folder:

```
$ ls -l ../bowtie2/*.bai
```

Create a directory for all of the output generated by deepTools, and move into that directory:

```
$ cd ~/ngs_course/chipseq/results/chip_qc		
 		
$ mkdir deeptools 		
 	
$ cd deeptools
```


Let's load the module and we are ready to get started:

```
$ module load seq/deeptools/2.2.0
```

### Calculation of the read coverage scores using the `multiBamSummary` tool

The `multiBamSummary` tool will calculate the read coverage scores for specific genomic regions between samples and provide the output as a binary compressed numpy array (.npz) file; however, the analysis can be performed on the entire genome by changing the mode of this tool to `bins`. It will also output a `readCounts.tab` file that contains a list read counts per sample for every 10,000bp region in the genome. 

```
$ multiBamSummary bins --ignoreDuplicates -p 6 \
--bamfiles ../../bowtie2/*aln.bam \
-out deeptools_multiBAM.out.npz \
--outRawCounts readCounts.tab
```

### Visualizing read coverage quality metrics

Now that we have the read coverage scores calculated for all samples, we can now analyze the coverage between samples using a variety of the *deepTools* tools:

#### 1. Sample correlation - `plotCorrelation` tool

The `plotCorrelation` tool allows us to visualize the similarity between samples based on their read coverage of regions of the genome. For example, we can compare two samples to determine whether they have similar coverage profiles with either a heatmap or a scatterplot:

![correlate](../img/QC_bamCorrelate_deeptools.png)

```
$ plotCorrelation --corData deeptools_multiBAM.out.npz \
--plotFile deepTools_scatterplot.png \
--corMethod pearson \
--whatToPlot scatterplot \
--labels Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2
```
<img src="../img/deepTools_scatterplot.png" width="600">

We expect high correlations between replicates, and lower correlations between samplegroups. However, we do not observe this when looking at read coverage on chromosome 12. Specifically, we see that Input-Rep1 does not correlate well with any of the other samples. If this were for the entire genome, we might be concerned that we would not have reporducibility between replicates for many of the peaks, and that Input-Rep1 is a potential outlier.

The same `plotCorrelation` tool can be used to examine the read coverage similarity using a heatmap to perform heirarchical clustering and determine whether our sample groups cluster well (i.e. have similar read coverage profiles within and between sample groups). The lack of correlation between replicates is even more visible in the heatmap, as is the lack of correlation with Input-Rep1 and all other samples.

```
$ plotCorrelation --corData deeptools_multiBAM.out.npz \
--plotFile deepTools_heatmap.png \
--corMethod pearson \
--whatToPlot heatmap \
--labels Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2 \
--plotNumbers
```
<img src="../img/deeptools_heatmap.png" width="600">

#### 2. Sample variability - `plotPCA` tool

The next quality metric we will explore is the principal component analysis (PCA) of our read coverage calculations. PCA can be used to determine whether samples display greater variability between experimental conditions than between replicates of the same treatment based on information (read coverage values) from thousands of regions. PCA is also useful to identify unexpected patterns, such as those caused by batch effects or outliers. 

You will use the tool `plotPCA` to sort the principal components according to the amount of variability of the data that they explain and generate two plots:

- the PCA plot for the top two principal components eigenvalues 
- the Scree plot for the top five principal components where the bars represent the amount of variability explained by the individual factors and the red line traces the amount of variability is explained by the individual components in a cumulative manner [[1]](http://deeptools.readthedocs.org/en/latest/content/tools/plotPCA.html)

![PCA](../img/PCA_deeptools.png)

```
$ plotPCA --corData deeptools_multiBAM.out.npz \
--plotFile deepTools_pcaplot.png \
-T "PCA of read counts" \
--outFileNameData deeptools_pcaProfile.tab \
--labels Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2
```

![pca_mov10](../img/deepTools_pcaplot.png)

Similar to the correlation plots, we see little clustering of the replicates. The variation between samplegroups does not account for the major sources of variation in the data. If this dataset were for the entire genome, I would be concerned that I would return few peaks that would be agreed upon by both replicates.

#### 3. Sample sequencing depth - `plotCoverage` tool

The `plotCoverage` tool will generate plots to explore the average number of reads per base pair in the genome. The tool will generate two plots, giving the frequencies of read coverage and the fraction of bases versus read coverage.

![coverage](../img/plotCoverage_deeptools.png)

```
$ plotCoverage --bamfiles ../../bowtie2/*aln.bam \
--region chr12:1,000,000-33,800,000 \
--ignoreDuplicates \
--plotFile deepTools_coverageplots.png \
--labels Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2
```

![coverage_mov10](../img/deepTools_coverageplots.png)

The mean coverage is very low, with only mean coverage per bp in region < 0.5 for all samples. Less than 10% of all positions have a read coverage of 2.0. If this were a full dataset, I would be worried that the coverage may be too low to accurately call peaks. However, we expect this low coverage with our subsetted data for chromosome 12.

#### 4. Sample signal strength - `plotFingerprints` tool

The `plotFingerprints` tool "determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample" [[2](http://deeptools.readthedocs.org/en/latest/content/tools/plotFingerprint.html)].  

"For factors that will enrich well-defined, rather narrow regions (e.g. transcription factors such as p300), the resulting plot can be used to assess the strength of a ChIP, but the broader the enrichments are to be expected, the less clear the plot will be" [[2](http://deeptools.readthedocs.org/en/latest/content/tools/plotFingerprint.html)].

The tool will generate a plot for the cumulative read coverages for each sample.

![fingerprint](../img/plotFingerprint_deeptools.png)

```
$ plotFingerprint \
--bamfiles ../../bowtie2/*aln.bam \
--minMappingQuality 25 --skipZeros \
--region chr12:1,000,000-33,800,000 \
--numberOfSamples 50000 \
-T "Fingerprints of different samples"  \
--plotFile deeptools_fingerprints.png \
--outRawCounts deeptools_fingerprints.tab
```

![fingerprint_mov10](../img/deepTools_fingerprints.png)

The input samples show the least enrichment, which is good, but Pou5f1 rep2 and Nanog rep1 exhibit little enrichment compared to input control. Also of note is that some of the samples don't have any coverage for ~20-50% of the region. While fingerprint plots are a nice way to look at signal strength, I would put more weight in the `phantompeakqualtools` results to determine whether we can sufficiently detect signal to noise in our samples. 

### Summary of ChIP-QC metrics

In summary, the dataset seems to have good signal-to-noise (based on phantompeakqualtools, though we would need to check values for all samples), indicating that we should be able to accurately call peaks in our dataset. That being said, the dataset exhibits extremely low coverage, indicating that if we had higher coverage, we would probably be able to call more peaks. In addition, the replicates do not cluster very well, so reproducibility between replicates may be an issue, which could result in few agreed upon peak calls.

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

