# Session IV: Single-cell RNA-seq Homework


## Running the single-cell RNA-seq workflow 

> **NOTE: When working on O2, be sure to run this in `/n/scratch2/` rather than your home directory .**

**Q1**.  Create a directory on `/n/scratch2` using your eCommons user ID as the directory name, if you have not already done so. Within that directory create a directory called `scRNA-seq`.

**Q2**. We are exploring public data from [Kang et al, 2017](https://www.nature.com/articles/nbt.4042), which was a study designed "to characterize the cell-type specificity and inter-individual variability of response to IFN-β, a potent cytokine that induces genome-scale changes in the transcriptional profiles of immune cells". From each of eight lupus patients, PBMCs were activated with recombinant IFN-β or left untreated for 6 h, a time point we previously found to maximize the expression of interferon-sensitive genes in dendritic cells and T cells. Two pools, IFN-β-treated and control, were prepared with the same number of cells from each individual and loaded onto the 10× Chromium instrument. A total of 14,619 (control) and 14,446 (stimulated) cell-containing droplets were obtained.

**We have downloaded this data** and made it available for you on O2.

**a.**  Setup a project directory structure within the `scRNA-se`q directory as shown below and copy over the count matrix for the stimulated condition into the appropriate directory from: `/n/groups/hbctraining/ngs-data-analysis-longcourse/sessionVI_hmwk/immune_stimulated_expression_matrix.txt.gz`

```bash
├── scRNA-seq
      ├── data/
      ├── logs/
      ├── meta/
      ├── scripts/
      ├── results/
      
```

> **Note**: We have both stimulated and control matrices available if you want to try to analyze more than a single condition.

**b.** Get O2 set-up to run the single-cell RNA-seq workflow

* Get X11 set-up on your laptop:
    * In class we ran single-cell RNA-seq on our laptops within RStudio after copying over the required files. If you are having trouble with any of the steps below with respect to X11, you can run it locally. 
    * The setup instructions are on [this page](https://wiki.rc.hms.harvard.edu/display/O2/Using+X11+Applications+Remotely#SSH+X11+Forwarding) under the "SSH X11 Forwarding" subheader.

> **Note 1**: Please note that once you set up the X11 forwarding, you will have use additional arguments when you use "ssh" to log on to O2 and  "srun" to start an interactive session. These are outlined clearly in the instructions linked above. 

> **Note 2**: Please test that it is working with a small example, maybe a simple ggplot2 code chunk or by using the example() function with another function that makes plots.

> **Note 3**: If you run into problems please [reach out to HMSRC folks](rchelp@hms.harvard.edu), this is known to be tricky! 

* Set up your personal libraries for R:
    * Check your .Renviron file (this should be in your home directory) and see what your R_LIBS_USER variable is set to. You will want it to be: `R_LIBS_USER="/n/groups/hbctraining/R/R-3.5.1/library"`
    * If something already exists in this file comment it out with a # sign.

> **Note**: If you didn't attend Session VI  and don't have an .Renviron file please [see this lesson](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/R_automation.html).

**c.** Run the single-cell RNA-seq workflow

* Perform quality control analysis on the counts. To do so, you will first need to read in the counts with the code listed below.

* Second, you will need to compute various QC metrics as we did in class (nUMI, nGene, log10GenesPerUMI, etc)

> **Note**: The mitochondrial genes appear to be absent in this dataset, so do not worry about calculating or plotting the mitochondrial ratio metrics.


```r

# Read in data
stim.data <- read.table("~/data/immune_stimulated_expression_matrix.txt.gz", 
sep = "\t")

# Set up stimulated object
stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim@meta.data$stim <- "STIM"

# Extract metadata
metadata_hw <- stim@meta.data

# Extract raw counts
counts_hw <- stim@raw.data

# Create a sparse matrix for more efficient computation
counts_hw <- as(as.matrix(counts_hw), "dgCMatrix")

```

*What thresholds did you use to filter the cells?*

* Determine the variable genes for your data, adjust parameters as needed.
    * *What number of variable genes were returned?*
    * *What was the full function that you used?*

* Determine whether to regress out the cell cycle phase, then regress out the uninteresting sources of variation. Determine the number of PCs to use for downstream clustering.
    * *What number of PCs will you use to determine the cell clusters?*
    * *Upload the PCElbowPlot.*

* Perform the clustering - use the cell type markers below to explore the presence of expected cell types and the resolution for the clustering:
    * **CD4 T cells**: IL7R, CD3D
    * **CD4 Naive T cells**: SELL
    * **CD4 Memory T cell**s: CREM
    * **CD14+ Monocytes**: CD14, LYZ, CCL2
    * **CD8 T cells**: CD8A
    * **B cell**s: MS4A1, CD79A (unactivated)
    * **FCGR3A+ Monocytes**: FCGR3A, MS4A7
    * **NK cells**: GNLY
    * **Dendritic cells**: FCER1A, CST3 
    * **Megakaryocytes**: PPBP

*What was the final resolution chosen for the clustering?*

* Perform marker identification; use the markers and previous known marker types to identify the cell type of each cluster.
    * *Upload a TSNE plot with the clusters labelled by cell type.*
