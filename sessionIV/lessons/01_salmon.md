---
title: "Quantitation of transcript abundance using Salmon"
author: "Mary Piper and Meeta Mistry"
date: "Thursday, January 19th, 2017"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 1.25 hours

## Learning Objectives

* Explore using lightweight algorithms to quantify reads to pseudocounts
* Understand how Salmon performs quasi-mapping and transcript abundance estimation
* Learn how to perform differential gene expression on pseudocounts

## Alignment-free quantification of gene expression

In the standard RNA-seq pipeline that we have presented so far in this course, we have taken our reads post-QC and aligned them to the genome using our transcriptome annotation (GTF) as guidance. The goal is to identify the genomic location where these reads originated from. Another strategy for quantification which has more recently been introduced involves **transcriptome mapping**. Tools that fall in this category include [Kallisto](https://pachterlab.github.io/kallisto/about), [Sailfish](http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html) and [Salmon](https://combine-lab.github.io/salmon/); each working slightly different from one another. (For this course we will explore Salmon in more detail.) Common to all of these tools is that we **avoid base-to-base alignment of the reads**, which is a time-consuming step, and these tools **provide quantification estimates much faster than do standard approaches** (typically 20 times faster) with improvements in accuracy. These estimates, often referred to as 'pseudocounts' are then converted for use with DGE tools like DESeq2. 

<img src="../img/alignmentfree_workflow_june2017.png" width="500">

### What is Salmon?

[Salmon](http://salmon.readthedocs.io/en/latest/salmon.html#using-salmon) is based on the philosophy of lightweight algorithms, which uses the sequence of genes or transcripts as input (in FASTA format) and do not align the whole read. However, many of these lightweight tools, in addition to the more traditional alignment-based tools, lack sample-specific bias models for transcriptome-wide abundance estimation. Sample-specific bias models are needed to account for known biases present in RNA-Seq data including:

- GC bias
- positional coverage biases
- sequence biases at 5' and 3' ends of the fragments
- fragment length distribution
- strand-specific methods

If not accounted for, these biases can lead to unacceptable false positive rates in differential expression studies [[1](http://salmon.readthedocs.io/en/latest/salmon.html#quasi-mapping-based-mode-including-lightweight-alignment)]. The **Salmon algorithm learns these sample-specific biases and accounts for them in the transcript abundance estimates** by combining a quasi-mapping approach with a dual inference method that is used to account for sample-specific biases and parameters. Salmon is extremely fast at "mapping" reads to the transcriptome and often more accurate than standard aproaches [[2](http://salmon.readthedocs.io/en/latest/salmon.html#quasi-mapping-based-mode-including-lightweight-alignment)]. 

### How does Salmon estimate transcript abundances?
Similar to standard base-to-base alignment approaches, the quasi-mapping approach utilized by Salmon requires a reference index to determine the position and orientation information for where the fragments best "map" prior to quantification [[1](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985/RapMap-a-rapid-sensitive-and-accurate-tool-for)]. 

#### **Indexing** 

This step involves creating an index to evaluate the sequences for all possible unique sequences of length k (kmer) in the **transcriptome** (genes/transcripts).

**The index helps creates a signature for each transcript in our reference transcriptome.** The Salmon index has two components:

- a suffix array (SA) of the reference transcriptome
- a hash table to map each transcript in the reference transcriptome to it's location in the SA (is not required, but improves the speed of "mapping" drastically)

#### **Quasi-mapping and quantification** 

The quasi-mapping approach estimates the numbers of reads mapping to each transcript, then generates final transcript abundance estimates after modeling sample-specific parameters and biases. The quasi-mapping approach is described below, with details provided by the Rapmap tool [[3,4]](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985/RapMap-a-rapid-sensitive-and-accurate-tool-for), which provides the underlying algorithm for the quasi-mapping.

- **Step 1:** Determine best mapping for each read/fragment and estimate number of reads/fragments mapping to each transcript

	<img src="../img/salmon_quasialignment.png" width="750">
	
	
	>RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes. A. Srivastava, H. Sarkar, N. Gupta, R. Patro. Bioinformatics (2016) 32 (12): i192-i200.
	
	As detailed in the figure above, the quasi-mapping procedure performs the following steps [[1](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985/RapMap-a-rapid-sensitive-and-accurate-tool-for)]:

	1. The read is scanned from left to right until a k-mer that appears in the hash table is discovered.
	2. The k-mer is looked up in the hash table and the SA intervals are retrieved, giving all suffixes containing that k-mer
	3. The maximal matching prefix (MMP) is identified by finding the longest read sequence that exactly matches the reference suffixes.
	4. We could search for the next MMP at the position following the MMP, but often natural variation or a sequencing error in the read is the cause to the mismatch from the reference, so the beginning the search at this position would likely return the same set of transcripts. Therefore, Salmon identifies the next informative position (NIP), which is the next position in the read where the SA search is likely to return a different set of transcripts than those returned for the previous MMP.
	5. This process is repeated until the end of the read.
	6. The final mappings are generated by determining the set of transcripts appearing in all MMPs for the read. The transcripts, orientation and transcript location are output for each read.
	
>
> **NOTE:** that if there are k-mers in the reads that are not in the index they are not counted. As such, trimming is not required when using this method.
	
	
- **Step 2:** Adjust abundance estimates based on RNA-Seq biases and sample-specific parameters

	Using the estimates of transcript abundance generated during the quasi-mapping procedure, a dual inference method is used to account for sample-specific biases and parameters. The biases are learned by one phase of the algorithm, which models the fragment-transcript assignment scores based on the following [[5](http://biorxiv.org/content/biorxiv/early/2016/08/30/021592.full.pdf)]:

	- chance of observing a fragment length given a particular transcript (derived from fragment length distribution)
	- chance fragment starts at particular position on transcript
	- concordance of fragment orientation based on specific protocol (stranded, paired-end)
	- chance fragment came from transcript based on score obtained from BAM (if BAMs were used instead of FASTQ)
	- positional, sequence-specific, and GC content based on computed mappings

	The model continuously learns and updates the transcript abundance estimates, and the the second phase of the algorithm refines the estimates by using the expectation maximization (EM) or variational Bayes optimization (VBEM)  [[5](http://biorxiv.org/content/biorxiv/early/2016/08/30/021592.full.pdf)]. The maximum likelihood estimates output from the EM or VBEM factorized likelihood function represent the estimated number of fragments derived from each transcript.
	
> **NOTE:** To have Salmon correct for these biases you will need to specify the appropriate parameters when you run it. Before using these parameters it is advisable to assess your data using tools like [Qualimap](http://qualimap.bioinfo.cipf.es/) to look specifically for the presence of these biases in your data and decide on which parameters would be appropriate. 
> 
> To correct for the various sample-specific biases you could add the following parameters to the sailfisg command:
>
> * `--seqBias` will enable it to learn and correct for sequence-specific biases in the input data.
> * `--gcBias` to learn and correct for fragment-level GC biases in the input data
> * `--posBias` will enable modeling of a position-specific fragment start distribution
>


## Running Salmon on Orchestra

First start an interactive session and create a new directory for our Salmon analysis:

```bash
$ bsub -Is -q interactive bash

$ mkdir ~/ngs_course/rnaseq/salmon

$ cd ~/ngs_course/rnaseq/salmon
```   
> Salmon is not available as a module on Orchestra, but it is installed as part of the bcbio pipeline. Since we already have the appropriate path (`/opt/bcbio/centos/bin`) in our `$PATH` variable we can use it by simply typing in `salmon`.     

As you can imagine from the description above, when running Salmon there are also two steps.

**Step 1: Indexing**
 "Index" the transcriptome (transcripts or genes) using the `index` command:
    
```bash
## DO NOT RUN THIS CODE
$ salmon index -t transcripts.fa -i transcripts_index --type quasi -k 31
```
> **NOTE:** Default for salmon is --type quasi and -k 31, so we do not need to include these parameters in the index command. The kmer default of 31 is optimized for 75bp or longer reads, so if your reads are shorter, you may want a smaller kmer to use with shorter reads (kmer size needs to be an odd number).
> 
**We are not going to run this in class, but it only takes a few minutes.** We will be using an index we have generated from transcript sequences (all known transcripts/ splice isoforms with multiples for some genes) for human. The transcriptome data (FASTA) was obtained from the [Ensembl ftp site](ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz).


**Step 2: Quantification:**
Get the transcript abundance estimates using the `quant` command and the parameters described below (more information on parameters can be found [here](http://salmon.readthedocs.io/en/latest/salmon.html#id5)):

* `-i`: specify the location of the index directory; for us it is `/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl37.idx/`
* `-l SR`: library type - specify stranded single-end reads (more info available [here](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype))
* `-r`: list of files for sample
* `--useVBOpt`: use variational Bayesian EM algorithm rather than the ‘standard EM’ to optimize abundance estimates (more accurate) 
* `-o`: output quantification file name
* `--writeMappings=salmon.out`: instead of printing to screen write to a file
   
To run the quantification step on a single sample we have the command provided below. Let's try running it on our subset sample for `Mov10_oe_1.subset.fq`:

```bash
$ salmon quant -i /groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl37.idx/ \
 -l SR \
 -r ~/ngs_course/rnaseq/raw_data/Mov10_oe_1.subset.fq \
 -o Mov10_oe_1.subset.salmon \
 --writeMappings=salmon.out \
 --useVBOpt 
```

>**NOTE:** Paired-end reads require both sets of reads to be given in addition to a [paired-end specific library type](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype):
`salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant`

## Salmon output

You should see a new directory has been created that is named by the string value you provided in the `-o` command. Take a look at what is contained in this directory:

    $ ls -l Mov10_oe_1.subset.salmon/
    
There is a logs directory, which contains all of the text that was printed to screen as Sailfish was running. Additionally, there is a file called `quant.sf`. 

This is the **quantification file** in which each row corresponds to a transcript, listed by Ensembl ID, and the columns correspond to metrics for each transcript:

```bash
Name    Length  EffectiveLength TPM     NumReads
ENST00000632684.1       12      3.00168 0       0
ENST00000434970.2       9       2.81792 0       0
ENST00000448914.1       13      3.04008 0       0
ENST00000415118.1       8       2.72193 0       0
ENST00000631435.1       12      3.00168 0       0
ENST00000390567.1       20      3.18453 0       0
ENST00000439842.1       11      2.95387 0       0

....

```

*  The first two columns are self-explanatory, the **name** of the transcript and the **length of the transcript** in base pairs (bp). 
*  The **effective length** represents the the various factors that effect the length of transcript due to technical limitations of the sequencing platform.
* Salmon outputs ‘pseudocounts’ which predict the relative abundance of different isoforms in the form of three possible metrics (KPKM, RPKM, and TPM). **TPM (transcripts per million)** is a commonly used normalization method as described in [[1]](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2820677/) and is computed based on the effective length of the transcript.
* Estimated **number of reads** (an estimate of the number of reads drawn from this transcript given the transcript’s relative abundance and length)

 
## Running Salmon on multiple samples 

We just ran Salmon on a single sample (and keep in mind a subset of chr1 from the original data). To obtain meaningful results we need to run this on **all samples for the full dataset**. To do so, we will create a shell script which will submit each Salmon run as a job to Orchestra.

Open up a script in `vim`:

	$ vim salmon_all_samples.sh

Now we can create a for loop to iterate over all FASTQ samples, and submit a job to **run Salmon on each sample in parallel**. We begin by listing all BSUB directives to specify the resources we are requesting including memory, cores and wall time.

Next comes the Salmon command. Note, that we are adding a parameter called `--numBootstraps` to the Salmon command. Salmon has the ability to optionally compute bootstrapped abundance estimates. **Bootstraps are required for estimation of technical variance**. Bootstrapping essentially takes a different sub-sample of reads for each bootstapping run for estimating the transcript abundances. The technical variance is the variation in transcript abundance estimates calculated for each of the different sub-samplings (or bootstraps). We will discuss this in more detail in the next lesson.

> *NOTE:* We are iterating over FASTQ files in the full dataset directory, located at `/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/full_dataset/`

> *NOTE:* We are using a new command called `basename` in our for loop. This command allows us to string split the file path to obtain a short name for filenaming purposes. We will discuss this in more detail when we talk about [Unix automation](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/blob/may2017/sessionIV/lessons/04_automating_workflow.md).


```bash
for fq in /groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/full_dataset/*.fastq
 do 
   base=`basename $fq .fastq`
   bsub -q mcore -n 6 -W 1:30 -R "rusage[mem=4000]" -J $base.mov10_salmon -o %J.$base.out -e %J.$base.err \
   salmon quant -i /groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl37.idx/ \
   -p 6 -l SR -r $fq --useVBOpt --numBootstraps 30 -o $base.salmon
 done
```

Save and close the script. This is now ready to run. **We are not going to run this script in class**, since it might take awhile and will interfere with the files we have already generated for you for use with the statistical analysis below.

> *NOTE:* PC users will want to add the `--auxDir` parameter to the Salmon command and provide an alternate name for the directory. By default it will be named `aux` which interferes with the decompressing process when bringing files over locally to run `tximport/DESeq2` in R.  

## Performing DE analysis on Pseudocounts

The pseudocounts generated by Salmon are represented as normalized TPM (transcripts per million) counts and map to transcripts or genes, depending on the input in the index step; in our case it was transcripts. These **need to be converted into non-normalized count estimates for performing DESeq2 analysis**. In addition, we need to be able to use the resulting non-normalized values which contain decimal places as input to create the DESeq object. Finally, to use DESeq2 we need to be able to tell which transcript is associated with which gene, since DESeq2 performs gene-level differential expression.

### Setting up to run DESeq2 on pseudocount data:

You can download the directory with the quant.sf files for the 8 full datasets using the link below. Once you have them downloaded continue to follow the set of instructions:

1. [Download Salmon files](https://www.dropbox.com/s/aw170f8zge01jpq/salmon.zip?dl=0)
2. Decompress (unzip) the zip archive and move the folder to an appropriate location (i.e `~/Desktop`)
3. Open RStudio and select 'File' -> 'New Project'  -> 'Existing Directory' and navigate to the `salmon` directory 
4. Open up a new R script ('File' -> 'New File' -> 'Rscript'), and save it as `salmon_de.R`

Your Rstudio interface should look something like the screenshot below:

<img src="../img/salmon_rstudio.png" size=600>

The developers of DESeq2 have developed a package that can make the conversion of alignment-free methods of quantification compatible for DESeq2. This package is called [`tximport`](https://bioconductor.org/packages/release/bioc/html/tximport.html) and is available through Bioconductor. `tximport` imports transcript-level abundance, estimated counts and transcript lengths, and summarizes this into matrices for use with downstream gene-level analysis packages. 

**Step 1:** Install the `tximport` package and the `readr` package (you'll only need to do this once):
    
```R
# Install from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("readr")
```

**Step 2:** Load the required libraries:

```R
# Load libraries
library(tximport)
library(readr)
library(DESeq2)
library(biomaRt) # tximport requires gene symbols as row names

```

**Step 3:** Load the quantification data that was output from Salmon:

```R
## List all directories containing data  
samples <- list.files(path = ".", full.names = F, pattern="\\.salmon$")
    
## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")
    
## Since all quant files have the same name it is useful to have names for each element
names(files) <-  samples
 ```
    
> **OPTION 2: An alternative to this is having absolute paths instead of relative paths.** This would be useful so you can run this from anywhere in your filesystem.
>

```R
## DO NOT RUN
dir <- getwd()
files <- file.path(dir, samples, "quant.sf")
	
## Create your own function
assignNames <- function(x){
		strsplit(x, "/")[[1]][6]
		}
names(files) <- sapply(files, assignNames, USE.NAMES=F)
```

Either of these methods will work, or even a combination of the two. The **main objective here is to add names to our quant files which will allow us to easily discriminate between samples in the final output matrix**. 

**Step 4.** Create a dataframe containing Ensembl Transcript IDs and Gene symbols

Our Salmon index was generated with transcript sequences listed by Ensembl IDs, but `tximport` needs to know **which genes these transcripts came from**, so we need to use the `biomaRt` package to extract this information. However, since BiomaRt has been a little unreliable we are actually not going to use this code right now.

> *NOTE:* Keep in mind that the Ensembl IDs listed in our Salmon output contained version numbers (i.e ENST00000632684.1). If we query Biomart with those IDs it will not return anything. Therefore, before querying Biomart in R do not forget to strip the version numbers from the Ensembl IDs.


```R
## DO NOT RUN

# Create a character vector of Ensembl IDs		
ids <- read.delim(files[1], sep="\t", header=T)    # extract the transcript ids from one of the files
ids <- as.character(ids[,1])
require(stringr)
ids.strip <- str_replace(ids, "([.][0-9])", "")

# Create a mart object
# Note that we are using an archived host, since "www.ensembl.org" gave us an error
mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host="mar2016.archive.ensembl.org"))

# Get official gene symbol and Ensembl gene IDs
tx2gene <- getBM(
    filters= "ensembl_transcript_id", 
     attributes= c("ensembl_transcript_id", "external_gene_name"),
     values= ids.strip,
     mart= mart)
     
```

**We have already run the above code for you and saved the output in a text file which is in the salmon directory.** Load it in using: 

```R
tx2gene <- read.delim("tx2gene.txt",sep="\t")
```
    
**Step 5:** Run tximport to summarize gene-level information    
```R
?tximport   # let's take a look at the arguments for the tximport function

txi <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene, reader=read_tsv, ignoreTxVersion=TRUE)
```
### Output from `tximport`

The `txi` object is a simple list with three matrices: abundance, counts, length. 
```R
attributes(txi)
```
A final element 'countsFromAbundance' carries through the character argument used in the tximport call. The length matrix contains the average transcript length for each gene which can be used as an offset for gene-level analysis. 


### Using DESeq2 for DE analysis with pseudocounts

```R    
library(DESeq2)   # load this if you have not loaded it earlier

## Create a sampletable/metadata

# Before we create this metadata object, let's see what the sample (column) order of the counts matrix is:
colnames(txi$counts)

condition=factor(c(rep("Ctl",3), rep("KD", 2), rep("OE", 3)))
sampleTable <- data.frame(condition, row.names = colnames(txi$counts))

## Create a DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
```

Now you have created a DESeq object to proceed with DE analysis as we discussed in the last session!

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*





