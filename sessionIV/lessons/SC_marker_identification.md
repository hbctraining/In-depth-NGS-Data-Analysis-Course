# Single-cell RNA-seq marker identification

Now that we have the single cells clustered based on different cell types,  we are ready to move forward with identifying cluster markers. Seurat has the functionality to perform a variety of analyses for marker identification; for instance, we can identify markers of each cluster relative to all other clusters by using the `FindAllMarkers()` function. This function essentially performs a differential expression test of the expression level in a single cluster versus the average expression in all other clusters.

To be identified as a cluster or cell type marker, we can specify thresholds for the minimum percentage of cells expressing the gene in either of the two groups of cells and minimum differential expression between the two groups. 

Using Seurat for marker identification is a rather quick and dirty way to identify markers. Usually the top markers are relatively trustworthy; however, because of inflated p-values, many of the less significant genes are not so trustworthy as markers. 

When looking at the output, we suggest looking for marker genes with large differences in expression between `pct.1` and `pct.2` and larger fold changes. For instance if `pct.1` = 0.90 and `pct.2` = 0.80, I might not be as excited about that marker. However, if `pct.2` = 0.1 instead, then I would be much more excited about it. Also, I look for the majority of cells expressing marker in my cluster of interest. If `pct.1` is low, such as 0.3, I again might not be as interested in it.

The results table contains the following columns:

- **p_val:** p-value not adjusted for multiple test correction
- **avg_logFC:** average log2 fold change. Positive values indicate that the gene is more highly expressed in the cluster.
- **pct.1**: The percentage of cells where the gene is detected in the cluster
- **pct.2**: The percentage of cells where the gene is detected on average in the other clusters
- **p_val_adj:** Adjusted p-value, based on bonferroni correction using all genes in the dataset, used to determine significance
- **cluster:** identity of cluster
- **gene:** Ensembl gene ID
- **symbol:** gene symbol
- **biotype:** type of gene
- **description:** gene description


```r
FindAllMarkers(seurat), 
```
