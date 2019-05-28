# Session V: ChIP-seq Homework

## Running the complete ChIP-seq workflow 
> **NOTE:** When working on O2, be sure to run this in `/n/scratch2/` rather than your home directory.

**Q1.**  Create a directory on `/n/scratch2` using your eCommons user ID as the directory name. Within that directory create a directory called `HCFC1_chipseq`.

**Q2.** You have strong evidence that HCFC1 is the transcription co-factor that associates with your protein of interest. To confirm this hypothesis you need to find binding regions for HCFC1 and see if they overlap with your current list of regions. The ENCODE project has ChIP-seq data for HCFC1 using a human liver cancer cell line HepG2 which contains 32 bp single end reads.  **We have downloaded this data and made it available for you on O2**. (NOTE: If you are interested in finding out more about the dataset you can find the [ENCODE record here](https://www.encodeproject.org/experiments/ENCSR529JYA/).

**a.**  Setup a project directory structure within the `HCFC1_chipseq` directory as shown below and copy over the raw FASTQ files from `/n/groups/hbctraining/ngs-data-analysis-longcourse/chipseq/HCFC1` into the appropriate directory: 

```bash
├── HCFC1_chipseq/
      ├── raw_data/
      ├── logs/
      ├── meta/
      ├── scripts/
      ├── results/
         ├── fastQC/
         ├── bowtie2/
         ├── unique_align/
         ├── macs2/
         └── sambamba/
         
```

**b.** Create a shell script that uses positional parameters and takes in a `.fastq` file as input and will perform the following:

* i. Run FASTQC 

* ii. Align reads with Bowtie2 using the parameters we used in class (Links to an external site.)Links to an external site..

> For the Bowtie2 index you will need to point to the hg19 index files from: `/n/groups/shared_databases/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/`

* iii. Change alignment file format from SAM to BAM (can be done using samtools or sambamba)

* iv. Sort the BAM file by read coordinate locations (can be done using sambamba or with samtools)

* v. Filter to keep only uniquely mapping reads (this will also remove any unmapped reads and duplicates) using sambamba

* vi. Index the final BAM file. This will be useful for visualization and QC.

> **Note 1: The script will require positional parameters, and using `basename` will help with the naming of output files**

```bash
#!/bin/bash
# Take the input FASTQ file and create a short name for output files
fq=$1
NAME=`basename $fq .fastq`

```

> **Note 2: make sure you load the required modules for each step**

> **Note 3: your final filtered BAM files should be located in `unique_align`**

**c.** Create a separate job submission script to **run the shell script you created in b) on all four `.fastq` files**. You have the option of running this in serial or parallel. Take a look at the [automation lesson from RNA-seq](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/06_automating_workflow.html#running-the-script-to-submit-jobs-in-parallel-to-the-slurm-scheduler) to help with setting up the job submission script.

**d.** Once the job submission script has finished running, use MACS2 to call peaks for each replicate using the input and treat BAM files and setting a q-value threshold of 0.05. *How many peaks are being called for each replicate?*

**e.** Run ChIPQC on O2. Create a QC HTML report. What can you conclude about the quality of the data based on these quality measures (individually, within group, and between group)?

> In class we ran ChIPQC on our laptops within RStudio after copying over the required files. If you are having trouble with any of the steps below with respect to X11, you can copy over the BAM files and peak calls to your laptop and run it locally. However, since this is a full dataset the BAM files are significantly larger which will mean that this will take significantly longer and you may run into problems due to lack of resources (disk space and/or memory).

In order to run ChIPQC on O2 you will first need to do the following:

* i. Get X11 setup for your laptop. The setup instructions are on [this page](https://wiki.rc.hms.harvard.edu/display/O2/Using+X11+Applications+Remotely#SSH+X11+Forwarding) under the "SSH X11 Forwarding" subheader.

> **Note 1:** Please note that once you set up the X11 forwarding, you will have use additional arguments when you use "ssh" to log on to O2 and  "srun" to start an interactive session. These are outlined clearly in the instructions linked above. 

> **Note 2:** Please test that it is working with a small example, maybe a simple `ggplot2` code chunk or by using the `example()` function with another function that makes plots.

> **Note 3:** If you run into problems please [reach out to HMSRC folks](rchelp@hms.harvard.edu), this is known to be tricky! 

* ii. Set up your personal libraries for R. If you attended SessionVI you will already have this done. If you didn't please see [this lesson](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/R_automation.html) and follow the instructions to setup a personal library.

* iii. Install the `ChIPQC` package into your personal library. **Remember that this is a Bioconductor package.**

* iv. Create a samplesheet for this dataset. A good idea would be to use the one from class as a template and modify accordingly. Use full paths to all of your files, just to be safe.

**f.** Sort each of the narrowPeak files using:

```bash
sort -k8,8nr HCFC1-rep1_peaks.narrowPeak > HCFC1-rep1_peaks_sorted.narrowPeak

sort -k8,8nr HCFC1-rep2_peaks.narrowPeak > HCFC1-rep2_peaks_sorted.narrowPeak
```

**g.** Use IDR to assess reproducibility between replicates using the sorted peaks in 2f above. Use the `awk` command from class to determine how many peaks meet an IDR cutoff of 0.05?

> **Note 1: You would normally set up the  IDR run on peaks called using a more liberal threshold, but for this question we will work with what we have.**

> **Note 2: Just perform the 1 step we did in class to generate the IDR stats; there is no need to do the remaining 2 steps in the IDR workflow.** (However, you can do them optionally if you are interested)

**h.**  These high confidence peaks from IDR can be used to explore downstream tools for functional enrichment and motif discovery. Use ChipSeeker to annotate the peaks and write those annotations to file. Take a look at binding site locations relative to the TSS, what does this tell us? *Hint: you can do this by looking at the summary, or you can use deepTools on the cluster to plot images that illustrate this.*

**i.** Use `clusterProfiler` to evaluate GO enrichment across the annotated peaks? What terms are over-represented?

**j.** **OPTIONAL:** Follow the links to lessons using the [MEME](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html#motif-analysis) suite of tools to explore motif analysis.
