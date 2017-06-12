---
title: "DEG analysis usin LRT in DESeq2"
author: "Meeta Mistry"
date: "October 17, 2016"
---

Approximate time: 35 minutes

## Learning Objectives 

* Introducing an alternative statistical test for differential expression analysis
* Extract results using the LRT and compare to Wald test
* Export results to file


## Hypothesis testing: Likelihood ratio test (LRT)

An alternative to pair-wise comparisons is to **analyze all levels of a factor at once**. By default the Wald test is used to generate the results table, but DESeq2 also offers the LRT which is used to identify any genes that show change in expression across the three levels. This type of test can be especially useful in analyzing time course experiments. 

To use the LRT, we use the `DESeq()` function but this time adding two arguments: 
1) to specify that we want to use the LRT `test`
2) the `reduced` model

```r
# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
```

Since our model only has one factor (`sampletype`), the reduced model is just the intercept. The LRT is comparing the full model to the reduced model to identify significant genes. **The p-values are determined solely by the difference in deviance between the full and reduced model formula (not log2 fold changes)**. 

Generally, this test will result in a larger number of genes than the individual pair-wise comparisons. While the LRT is a test of significance for differences of any level of the factor, one should not expect it to be exactly equal to the union of sets of genes using Wald tests (although we do expect a majority overlap).

Let's take a look at the results table:

```r
# Extract results
res_LRT <- results(dds_lrt)
```

You will find that similar columns are reported for the LRT test. One thing to note is, even though there are fold changes present they are not directly associated with the actual hypothesis test. Thus, when filtering significant genes from the LRT we use only the FDR as our threshold. *How many genes are significant at `padj < 0.05`?*

```r
sig_res_LRT <- subset(res_LRT, padj < padj.cutoff)
dim(sig_res_LRT)

# Get sig gene lists
sigLRT_genes <- rownames(sig_res_LRT)
length(sigLRT_genes)

# Compare with the sigOE and sigKD gene lists
sigOE_genes <- rownames(sigOE)
length(sigOE_genes)

sigKD_genes <- rownames(sigKD)
length(sigKD_genes)

sigOE_genes %in% sigLRT_genes
How many genes from the Mov10 overexpression Wald test are contained in the LRT gene set? And for the Mov10 knockdown? 

The number of significant genes observed from the LRT is quite high. We are **unable to set a fold change criteria here since the statistic is not generated from any one pairwise comparison.** This list includes genes that can be changing in any number of combinations across the three factor levels. It is advisable to instead increase the stringency on our criteria and lower the FDR threshold.

***

**Exercise**

1. Using a more stringent cutoff of `padj < 0.001`, count how many genes are significant using the LRT method.
2. Set the variables `OEgenes` and `KDgenes`to contain the genes that meet the  threshold `padj < 0.001`.
3. Find the overlapping number of genes between these gene sets and the genes from LRT at `padj < 0.0001`.

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*

