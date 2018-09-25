# Single-cell RNA-seq clustering analysis

Now that we have our high quality cells, we want to know the different cell types present within our population of cells. To do this we are going to perform a clustering analysis. The workflow for this analysis is adapted from the following sources:

- Satija Lab: [Seurat v2 Guided Clustering Tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html)
- Paul Hoffman: [Cell-Cycle Scoring and Regression](http://satijalab.org/seurat/cell_cycle_vignette.html)

To identify clusters, the following steps will be performed:

1. **Normalization and transformation** of the raw gene counts per cell to account for **differences in sequencing depth** per cell.
2. Identification of high variance genes.
3. **Regression of sources of unwanted variation** (e.g. number of UMIs per cell, mitochondrial transcript abundance, cell cycle phase).
4. **Identification of the primary sources of heterogeneity** using principal component (PC) analysis and heatmaps.
5. **Clustering cells** based on significant PCs (metagenes).

To perform this analysis, we will be mainly using functions available in the Seurat package. Therefore, we need to load the Seurat library in addition to the tidyverse library.

```r
library(Seurat)
library(tidyverse)
```

To perform the analysis, Seurat requires the data to be present as a `seurat` object. To create the `seurat` object, we will be extracting the **filtered counts** and **metadata** stored in our `se_c` SingleCellExperiment object created during quality control. 
To access the counts from our SingleCellExperiment, we can use the `counts()` function:

```r
seurat_raw <- CreateSeuratObject(raw.data = counts(se_c),
                                 meta.data = colData(se_c) %>% data.frame())
```

### Normalizing counts, finding variable genes, and scaling the data

The first step in the analysis is to normalize the raw counts to account for differences in sequencing depth per cell. The raw counts are normalized using global-scaling normalization with the `NormalizeData()` function, which performs the following:

1. normalizes the gene expression measurements for each cell by the total expression 
2. multiplies this by a scale factor (10,000 by default)
3. log-transforms the result

```r
# Normalize counts for total cell expression and take log value                            
pre_regressed_seurat <- seurat_raw %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000)  
```

Following normalization, we want to identify the most variable genes to use for downstream clustering analyses. The `FindVariableGenes()` function can be called, which performs the following calculations:

1. calculates the average expression and dispersion for each gene
2. places these genes into bins
3. calculates a z-score for dispersion within each bin

This helps control for the relationship between variability and average expression. 

```r
# Find variable genes based on the mean-dispersion relationship based on z-score for dispersion. 
pre_regressed_seurat <-  pre_regressed_seurat %>%
                          FindVariableGenes(
                            mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            do.plot = FALSE)
```

It's recommended to set parameters as to mark visual outliers on dispersion plot - default parameters are for ~2,000 variable genes. There are some additional arguments, such as `x.low.cutoff`, `x.high.cutoff`, `y.cutoff`, and `y.high.cutoff` that can be modified to change the number of variable genes identified. Generally, we might be a bit concerned if we are returning 500 or 4,000 variable genes.

We can check the number of variable genes to see if it meets expectations. 

```r
# Check number of variable genes to determine if correct parameters used  
length(x = pre_regressed_seurat@var.genes)
```

We can plot dispersion (a normalized measure of to cell-to-cell variation) as a function of average expression for each gene to **identify a set of high-variance genes**. To check that the dispersions behave as expected, decreasing with increasing mean, and to identify the most variable genes, we can visualize the dispersions with the `VariableGenePlot()` function.

```r
# Plot variable genes
VariableGenePlot(pre_regressed_seurat)
```

The identified variable genes are going to be the genes used to **identify significant principal components** used to determine the **how similar individual cells are to each other for clustering analysis**. 

However, to identify the significant principal components the expression values need to be centered and scaled. Centering each gene will center the expression of each gene by subtracting the average expression of the gene for each cell. Scaling will divide the centered gene expression levels by the standard deviation. To perform the centering and scaling, we can use Seurat's `ScaleData()` function.

```r
# Scale and center data
pre_regressed_seurat <- pre_regressed_seurat %>%
                        ScaleData(model.use = "linear")
```

### Examining sources of variation in the data

Your single-cell dataset likely contains "uninteresting" sources of variation. This can include technical noise, batch effects, and/or uncontrolled biological variation (e.g. cell cycle). Similar to bulk RNA-seq analysis, we can use PCA to identify these sources of variation, which can then be regressed out prior to further analysis.

### Cell cycle scoring

Cell cycle variation is a common source of uninteresting variation in single-cell RNA-seq data. To examine cell cycle variation in our data, we assign each cell a score, based on its expression of G2/M and S phase markers. 

> An overview of the cell cycle phases is given in the image below:
> 
> <img src="../img/cell_cycle.png" width="300">
> 	
> *Adapted from [Wikipedia](https://en.wikipedia.org/wiki/Cell_cycle) (Image License is [CC BY-SA 3.0](https://en.wikipedia.org/wiki/Wikipedia:Text_of_Creative_Commons_Attribution-ShareAlike_3.0_Unported_License))*
> 
> - **G0:** Quiescence or resting phase. The cell is not actively dividing, which is common for cells that are fully differentiated. Some types of cells enter G0 for long periods of time (many neuronal cells), while other cell types never enter G0 by continuously dividing (epithelial cells).
> - **G1:** Gap 1 phase represents the **beginning of interphase**. During G1 there is growth of the non-chromosomal components of the cells. From this phase, the cell may enter G0 or S phase.
> - **S:** Synthesis phase for the replication of the chromosomes (also part of interphase).
> - **G2:** Gap 2 phase represents the **end of interphase**, prior to entering the mitotic phase. During this phase th cell grows in preparation for mitosis and the spindle forms.
> - **M:** M phase is the nuclear division of the cell (consisting of prophase, metaphase, anaphase and telophase).
	

At the HBC core, we have accumulated a nice list of genes associated with particular cell cycle phases. We are going to download the list of cell cycle phase marker genes by **right-clicking** [here](https://github.com/hbc/tinyatlas/raw/master/cell_cycle/Homo_sapiens.csv) and saving to the `data` folder.

To save the genes in the G2M and S phases as character vectors, we can subset the data frame:

```r
# Read in cell cycle genes
cell_cycle <- read.csv("data/Homo_sapiens.csv")

# Extract the G2/M genes
g2m_genes <- dplyr::filter(cell_cycle, phase == "G2/M") %>%
  pull(geneID) %>%
  as.character() 
  
# Extract the S genes
s_genes <- dplyr::filter(cell_cycle, phase == "S") %>%
  pull(geneID) %>%
  as.character() 
```

Now to score each gene for cell cycle, we can use Seurat's `CellCycleScoring()` function. The function scores cells based on their expression of the G2M and S phase marker genes, which should be anticorrelated in their expression levels, and cells expressing neither are likely not cycling and in G0/G1 phase. The `CellCycleScoring()` function stores S and G2/M scores in `seurat@meta.data` in the `S.Score` and `G2M.Score` columns, along with the predicted classification of each cell in either G2M, S or G1 phase in the `Phase` column.

```r
# Perform cell cycle scoring
pre_regressed_seurat <- CellCycleScoring(
  pre_regressed_seurat,
  g2m.genes = g2m_genes,
  s.genes = s_genes)
```

To determine whether the cells group by cell cycle, we can perform PCA using the cell cycle genes. If the cells group by cell cycle in the PCA, then we would want to regress out cell cycle variation. By default in the  `RunPCA()` function, the most variable genes identified previously are used to determine the PCs, but we can select specific genes using the `pc.genes` argument. We will use only the cell cycle genes to determine whether the cells cluster by cell cycle phase. 

```r
# Perform PCA and color by cell cycle phase
pre_regressed_seurat = RunPCA(
  pre_regressed_seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

PCAPlot(pre_regressed_seurat, 
        group.by= "Phase")
```

In our data, the cells don't really cluster by cell cycle, so we do not need to include `S.Score` and `G2M.Score` as variables for regression.

Before moving on to regressing out variation due to uninteresting sources, let's save the pre-regressed Seurat object so that we can come back to it later if needed:

```r
# Save pre-regression Seurat object
saveRDS(pre_regressed_seurat, 
        file = "data/seurat_pre_regress.rds")
```


>**NOTE:** Often we only want to analyze a subset of samples, cells, or genes. To subset the Seurat object, the `SubsetData()` function can be easily used. For example, to only cluster cells using a single sample group:
>
>```r
> pre_regressed_seurat <- SubsetData(seurat_raw, 
>                                    cells.use = rownames(seurat_raw@meta.data[which(seurat_raw@meta.data$interestingGroups == "control"), ])
>```

## Apply regression variables

Regressing variation due to uninteresting sources can improve downstream identification of principal components and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on the variables to regress.

To regress out these variables of uninteresting variation, we will use the `vars.to.regress` argument in the `ScaleData()` function. 

We generally recommend minimizing the effects of variable read count depth (`nUMI`) and mitochondrial gene expression (`mitoRatio`) as a standard first-pass approach. However, if the differences in mitochondrial gene expression represent a biological phenomenon that may help to distinguish cell clusters, then we advise not passing in `mitoRatio`.

When regressing out the effects of cell-cycle variation, include `S.Score` and `G2M.Score` in the `vars.to.regress` argument. Cell-cycle regression is generally recommended but should be avoided for samples containing cells undergoing differentiation.

In our data, the cell cycle phase did not appear to be a large source of variation in the data, so we do not need to regress it out. Therefore, we will only regress out variation due to `nUMI` and `mitoRatio`.

```r
# Define variables in metadata to regress
vars_to_regress <- c("nUMI", "mitoRatio")

# Regress out the uninteresting sources of variation in the data
seurat <- ScaleData(pre_regressed_seurat, 
                    vars.to.regress = vars_to_regress)
```


## Linear dimensionality reduction

Next, we perform principal component analysis (PCA) on the scaled data with `RunPCA()`. `ProjectPCA()` scores each gene in the dataset (including genes not included in the PCA) based on their correlation with the calculated components.

```r
# Perform the scoring for all genes
seurat <- seurat %>%
  RunPCA(do.print = FALSE) %>%
  ProjectPCA(do.print = FALSE)
```

## Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a "metagene" that combines information across a correlated gene set. Often it is useful to explore the PCs:

```r
# Print out the top 5 most variant genes (up and down) for top 5 PCs
PrintPCA(object = seurat, 
         pcs.print = 1:5, 
         genes.print = 5, 
         use.full = FALSE)
```

Determining how many PCs to include downstream is therefore an important step. To accomplish this, we plot the standard deviation of each PC as an elbow plot with our `plotPCElbow()` function.

```r
# Create elbow plot
PCElbowPlot(seurat)
```

Based on this plot, we can eye the plot, and the elbow appears to be around PC 7 or 8. While this gives us a good idea of the number of PCs to include, a more quantitative approach may be a bit more reliable.

PC selection — identifying the true dimensionality of a dataset — is an important step for our clustering analysis, but can be challenging/uncertain. While there are a variety of ways to choose a threshold, we're going to calculate where the principal components start to elbow by taking the larger value of:

1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
2. The point where the percent change in variation between the consequtive PCs is less than 0.1%.

We will start by calculating the first metric:

```r
# Determine percent of variation associated with each PC
pct <- seurat@dr$pca@sdev / sum(seurat@dr$pca@sdev) * 100

# Calculate cumulative percents for each PC
cum <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]

co1
```
The first metric returns PC18 as the PC matching these requirements. Let's check the second metric, which identifies the PC where the percent change in variation between consequtive PCs is less than 0.1%:

```r
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.

co2
```

The second metric returned PC8. Now, to determine the selection of PCs, we will use the minimum of the two metrics:

```r
# Minimum of the two calculation
pcs <- min(co1, co2) # change to any other number

pcs
```

Based on these metrics, for the clustering of cells in Seurat we will use the first **eight PCs** to generate the clusters.

## Cluster the cells

We can now use these significant PCs to determine which cells exhibit similar expression patterns for clustering. To do this, Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. Details on this clustering methods are available in the Seurat paper.
We will use the `FindClusters()` function to perform the graph-based clustering. The `resolution` argument that sets the "granularity" of the downstream clustering, will need to be optimized to the experiment, with increased values leading to a greater number of clusters. We find that setting this parameter between `0.6`-`1.2` typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The cluster IDs are saved in the `seurat@ident` slot.

We provide a series of resolution options during clustering, which can be used downstream to choose the best resolution.


```r
# Find cell clusters
seurat <- FindClusters(
  seurat,
  dims.use = 1:pcs,
  force.recalc = TRUE,
  print.output = TRUE,
  resolution = c(0.6, 0.8, 1.0, 1.2),
  save.SNN = TRUE)
```

## Creating t-SNE plots

Seurat continues to use t-distributed stochastic neighbor embedding (t-SNE) as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on t-SNE components, cells within the graph-based clusters determined above should co-localize on the t-SNE plot. This is because the t-SNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the t-SNE, we suggest using the same PCs as input to the clustering analysis, although computing the t-SNE based on scaled gene expression is also supported using the `genes.use` argument. **Note that distance between clusters on the t-SNE plots does not represent degree of similarity between clusters.**

```r
# Choose a resolution
seurat <- SetAllIdent(object = seurat, id = "res.0.8")

# Run the TSNE and plot
seurat <- RunTSNE(
  seurat,
  dims.use = 1:pcs,
  do.fast = TRUE)
```

```r
# Plot the TSNE
TSNEPlot(object = seurat)
```

Once a resolution has been chosen, a useful feature in Seurat is the ability to recall the parameters that were used in the latest function calls for commonly used functions. For `FindClusters()`, the authors provide the function `PrintFindClustersParams()` to print a nicely formatted formatted summary of the parameters that were chosen.

```r
PrintFindClustersParams(seurat)
```

```r
# Save clustered cells
saveRDS(seurat, file = file.path(data_dir, "pbmcs_seurat_tsne.rds"))
```

# Evaluating clustering

In order to determine whether our clustering and resolution are appropriate for our experiment, it is helpful to explore a handful of markers for each of the major cell types that we expect to be in our data and see how they segregate.

The `FeaturePlot()` function from seurat makes it easy to visualize a handful of genes using the gene IDs stored in the Seurat object. For example if we were interested in exploring known immune cell markers, such as:

- T cells: "ENSG00000010610", "ENSG00000153563"
- Monocytes: "ENSG00000131495", "ENSG00000150337"
- B cells: "ENSG00000177455"
- Natural killer cells: "ENSG00000149294"
- DC: "ENSG00000140678"
- CD14: ENSG00000170458
- CD16a: ENSG00000203747


```r
FeaturePlot(object = seurat, 
            features.plot = c("ENSG00000010610", "ENSG00000153563", "ENSG00000131495", "ENSG00000150337", "ENSG00000177455", "ENSG00000149294", "ENSG00000140678", "ENSG00000203747", "ENSG00000170458"))
```

# QC

To determine whether our clusters might be due to artifacts such as cell cycle phase or mitochondrial expression, it can be 




> **NOTE:** Most single-cell RNA-seq datasets are too big to work with on a personal laptop, so you will need to use R on O2. To do this requires establishing a personal R library with the appropriate libraries installed. More information about setting up personal libraries [is available](https://wiki.rc.hms.harvard.edu/display/O2/Personal+R+Packages) from HMS RC. In addition to a personal R library, the analysis on O2 can be difficult if you cannot view the results. To view plots/images output on O2 requires X11 forwarding, and how to enable X11 configuration on your computer [is also detailed](https://wiki.rc.hms.harvard.edu/display/O2/Using+X11+Applications+Remotely) by HMS RC.
