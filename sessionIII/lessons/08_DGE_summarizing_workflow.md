---
title: "Summary of DGE workflow"
author: "Mary Piper"
date: "June 8, 2017"
---

Approximate time: 15 minutes

## Learning Objectives 

* Understand the commands needed to run a complete differential expression analysis

## Summary of differential expression analysis workflow

We have detailed the various steps in a differential expression analysis workflow, providing theory with example code. To provide a more succinct reference for the code needed to run a DGE analysis, we have summarized the steps in an analysis below:

1. Count normalization:
	
	```r
	# Check that the row names of the metadata equal the column names of the **raw counts** data
	all(colnames(data) == rownames(meta))
	
	# Create DESeq2Dataset object
	dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ condition)
	
	# Estimate size factors
	dds <- estimateSizeFactors(dds)
	
	# Output normalized counts
	normalized_counts <- counts(dds, normalized=TRUE)
	```
	
2. Exploratory data analysis (PCA & heirarchical clustering) - identifying outliers and sources of variation in the data:
	
	```r
	# Transform counts for data visualization
	rld <- rlog(dds, blind=TRUE)
	
	# Plot PCA 
	plotPCA(rld, intgroup="condition")
	
	# Extract the rlog matrix from the object
	rld_mat <- assay(rld)
	
	# Compute pairwise correlation values
	rld_cor <- cor(rld_mat)
	
	# Plot heatmap
	pheatmap(rld_cor)
	```
	
3. Run DESeq2:

	```r
	# Create DESeq2 dataset
	dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ condition)
	
	# Run DESeq2 differential expression analysis
	dds <- DESeq(dds)
	```
	
4. Check the fit of the dispersion estimates:
	
	```r
	# Plot dispersion estimates
	plotDispEsts(dds)
	``` 

5. Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:

	```r
	res <- results(dds, contrast = c("sample_group", "level_to_compare", "base_level"))
	```

6. Output significant results:

	```r
	sig_res <- subset(as.data.frame(res), res_tableOE$padj < padj.cutoff & abs(res_tableOE$log2FoldChange) > lfc.cutoff)
	```

7. Visualize results (volcano plots, heatmaps, etc.)
