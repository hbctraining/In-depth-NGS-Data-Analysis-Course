# Quality Control Analysis

### Dataset

The dataset we will be working with is comprised of 2,700  Peripheral Blood Mononuclear Cells (PBMC) sequenced on  the Illumina NextSeq 500. This dataset is freely available from 10X Genomics and is used as part of the [Seurat tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html).


### Setting up the R environment

Create a new R project entitled `single_cell_rnaseq`. Then, create the following directories:

- data
- results
- figures

**Right-click** the link [here](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and download the data into the `data` folder.

We will need to navigate to the `data` folder and **click on the file `pbmc3k_filtered_gene_bc_matrices.tar.gz`** to decompress it. 

Explore the files created by clicking on the `filtered_gene_bc_matrices/hg19` folder:

- **`barcodes.tsv`:** cellular barcodes
- **`genes.tsv`:** list of genes
- **`matrix.mtx`:** counts assigned to each gene for each single cell


Finally, create an Rscript and type the following note:

```r
# Single-cell RNA-seq analysis with Seurat - QC
```

Save the Rscript as `quality_control.R`. Now, we can load the necessary libraries:

# Explain the input for QC: tagcounts, rownames, colnames files

# Define computationally efficient object dgcmatrix

# Create count matrix

```r
library(SingleCellExperiment)
library(Matrix)
library(tidyverse)

counts = readMM("~/orch/hbc/PIs/john_gaber/indrop_scrnaseq_patients/sc-human/final/2018-07-19_sc-human/tagcounts.mtx")
rownames = read.csv("~/orch/hbc/PIs/john_gaber/indrop_scrnaseq_patients/sc-human/final/2018-07-19_sc-human/tagcounts.mtx.rownames", header = F)[["V1"]] %>% as.character()
colnames = read.csv("~/orch/hbc/PIs/john_gaber/indrop_scrnaseq_patients/sc-human/final/2018-07-19_sc-human/tagcounts.mtx.colnames", header = F)[["V1"]] %>% 
    as.character() %>% 
    make.names()
counts =  as(counts, "dgCMatrix")
rownames(counts) = rownames
colnames(counts) = colnames

```

# Create metadata - switch [[]] to dollar sign

```r
metadata = data.frame(row.names=colnames, samples = colnames, stringsAsFactors = F)

metadata[["lane"]] = gsub("\\..+", "", metadata[["samples"]])
metadata[["nUMI"]] = colSums(counts)
metadata[["nGenes"]] = colSums(counts>0)
metadata[["log10GenesPerUMI"]] = log10(metadata$nGene) / log10(metadata$nUMI)
```
# Annotation of the mitochondrial genes - table with only one Ensembl ID per row (could just deterine Ensembl IDs for mito genes)

```r
# annotation was download from ensembl biomart to match the version GRCh38.92
# AnnotationHub can be used.
rows = read_tsv("annotation.tsv") %>% 
    janitor::clean_names() %>% 
    as.data.frame() %>% 
    dplyr::select(gene_id = gene_stable_id,
                  gene_name,
                  gene_description,
                  biotype = gene_type,
                  entrezid = ncbi_gene_id,
                  chrom = chromosome_scaffold_name) %>% 
    group_by(gene_id, biotype, chrom, gene_description) %>% 
    summarise(gene_name = paste(unique(gene_name), collapse = ","),
              entrezid = paste(unique(entrezid), collapse = ",")) %>% 
    mutate(gene_name=ifelse(gene_name=="", gene_id, gene_name)) %>% 
    as.data.frame()
# mit
rrna = rows %>% dplyr::filter(grepl("rRNA",biotype)) %>% .[["gene_id"]] %>% 
    intersect(., rownames) # this is empty, not used
trna = rows %>% dplyr::filter(grepl("tRNA",biotype)) %>% .[["gene_id"]] %>% 
    intersect(., rownames) # this is empty, not used
mt = rows %>% dplyr::filter(chrom == "MT") %>% .[["gene_id"]] %>% intersect(., rownames)

metadata[["mtUMI"]] = colSums(counts[mt,], na.rm = T)
metadata[["mtUMI"]][is.na(metadata[["mtUMI"]])] = 0
metadata[["mitoRatio"]] = metadata$mtUMI/metadata$nUMI

idx = which(metadata$nUMI>100)
counts_c = counts[, idx]
metadata_c = metadata[idx,]
```

# Turn into single cell experiment 

```r
se = SingleCellExperiment(assays=list(raw=counts_c), colData = metadata_c)
saveRDS(se, "data/se.rds")
```

## Read in data

```r
library(cowplot)
library(ggplot2)
library(scales)
library(tidyverse)
library(SingleCellExperiment)
# Load bcbioSingleCell object
se <- readRDS(params$bcb_file)
metrics = colData(se) %>% as.data.frame
```

## Reads per cell

### Proportional histogram

With the proportional histogram you hope to see all of the samples with peaks in relatively the same location between 10,000 and 100,000 reads per cell. A shoulder would be indicative of many poor quality cells. 
- can't do without the number of reads

## Cell counts

### Bar plot

The cell counts are determined by the number of unique cellular barcodes detected. 

You expect the number of unique cellular barcodes to be around the number of sequenced cells (determined in step 1) or greater due to some hydrogels having more than one cellular barcode. 

>During the **inDrop** protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.

```r
metrics %>% 
  ggplot(aes(x=lane)) + geom_bar() + ggtitle("NCells")
metrics %>% 
  ggplot(aes(x=lane, y=log10(nGenes))) + geom_boxplot() + ggtitle("NCells vs NGenes")

# change lane to sample
```

## UMI counts (transcripts) per cell

### Raw ridgeline

The UMI counts per cell should be generally above 500, although usable, it's still low if between 500-1000 counts. If UMIs per cell is 500-1000 counts, then the cells probably should have been sequenced more deeply. 

```r
metrics %>% 
    ggplot(aes(color=lane, x=nUMI)) + geom_density() + scale_x_log10() + geom_vline(xintercept = 500)
```
    
## Genes detected per cell

### Raw ridgeline

Seeing gene detection in the range of 500-5000 is normal for **inDrop** analysis. Similar expectations for gene detection as for UMI detection.

```r
metrics %>% 
    ggplot(aes(color=lane, x=nGenes)) + geom_density() + scale_x_log10() + geom_vline(xintercept = 200)
```

## UMIs vs. genes detected

### Line plot

Poor quality cells are likely to have low genes and UMIs per cell. Therefore, a poor sample is likely to have cells in the lower left of the graph. Good cells should exhibit both higher number of genes per cell and higher numbers of UMIs. We also expect similar lines with similar slopes for all samples.

```r
metrics %>% 
    ggplot(aes(x=nUMI, y=nGenes, color=mitoRatio)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_vline(xintercept = 800)+
      facet_wrap(~lane)
```

## Mitochondrial counts ratio

### Raw ridgeline

This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. Poor quality samples for mitochondrial counts would have larger peaks above the 0.1 mitochondrial ratio mark, unless it is expected based on sample type.

```r
metrics %>% 
    ggplot(aes(color=lane, x=mitoRatio)) + geom_density() + scale_x_log10() + geom_vline(xintercept = params$max_mito_ratio)
```

## Novelty

### Raw ridgeline

We can see the samples where we sequenced each cell less have a higher overall novelty, that is because we have not started saturated the sequencing for any given gene for these samples. Outlier cells in these samples might be cells that we have a less complex RNA species than other cells. Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric.

```r
metrics %>%
    ggplot(aes(x=log10GenesPerUMI, color = lane)) +
    geom_density()
```

# Filtering

```r
keep = metrics %>% dplyr::filter(nUMI > 500 , nReads > 3000) %>% .[["samples"]]
se_c = se[,keep]
metrics_clean = colData(se_c) %>% as.data.frame()
saveRDS(se_c, file = "data/se_filtered.rds")
```

# Repeat QC plots