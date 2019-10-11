# Session III: Advanced R and Differential Gene Expression Homework 


## Using DESeq2 for gene-level differential expression analysis 

**Q1**.  The metadata below describes an experiment that you have setup for RNA-seq analysis. Assume the situation in which all required libraries have been loaded and the metadata table below and associated count matrix have been loaded into R as `txi` and `meta`. Use the information in the table to answer the following questions.  

**NOTE: This is an exercise in thinking about running DESeq2. You do not need to run any code in R/RStudio. Refer to the materials/lessons from class to answer the following questions.**

| Genotype | Celltype | Batch |
|:-----------:|:----------:| :----------:| 
|sample1 |	Wt	| typeA	| second |
|sample2 |	Wt	| typeA	| second |
|sample3 |	Wt	| typeA	| first |
|sample4 |	KO	| typeA	| first |
|sample5 |	KO	| typeA	| first |
|sample6 |	KO	| typeA	| second |
|sample7 |	Wt	| typeB	| second |
|sample8 |	Wt	| typeB	| first |
|sample9 |	Wt	| typeB	| second |
|sample10 |	KO	| typeB	| first |
|sample11 |	KO	| typeB	| first |
|sample12 |	KO	| typeB	| second |


a. Provide the line of code used to create a `DESeqDataSet` object called `dds` in which “Genotype” is the factor of interest and "Celltype" and "Batch" are other contributing sources of variation in your data.

b. Provide the line of code required to run DESeq2 on `dds`.

c. Provide the line of code to create a dispersion plot.

d. Provide the line of code to perform a Wald test comparison of the Celltype categories `typeA` versus `typeB` (i.e the **fold changes reported should reflect gene expression changes relative to `typeB`**).

e. Provide the line of code required to shrink the fold changes.

f.  Provide the line of code to write the results of the Wald test above to a file called "results_celltype.csv".

g. Provide the line of code to subset the results to return those genes with adjusted p-value < 0.05 and an absolute logFold2Change > 1.

 
## Working with the DESeq2 results table 

**Q2**.  Using the `de_script.R` that we created in class for the differential expression analysis, change the thresholds for adjusted p-value and log fold change to the following values:

* padj.cutoff: 0.01

* lfc.cutoff: 1.5

Using these new cutoffs, perform the following steps:

a. Subset `res_tableOE` to only return those rows that meet the criteria we specified above (adjusted p-values < 0.01 and log fold changes >1.5). Save the subsetted table to a data frame called `sig_table_hw_oe`. **Write the code below**:

b. There is a DESeq2 function that summarizes how many genes are up- and down-regulated using our criteria for alpha=0.01. Use this on the `res_tableOE`. **Write the code you would use, and also, list how many genes are up- and down-regulated**.

c. Get the gene names from `sig_table_hw_o`e and save them to a vector called `sigOE_hw`. **Write the code below**:

d. Write the `sigOE_hw` vector of gene names to a file called `sigOE_hw.txt` using the write() function. Ensure the genes are listed in a single column. **Write the code below and upload `sigOE_hw.txt` file**.

e. Use the newly created genes list (`sigOE_hw`) to perform over-representation analysis using clusterProfiler. If you find enriched GO terms, upload the dotplot and the cnetplot (genes colored by fold change) write the code you would use to generate a dot plot and cnet plot .

 

## Visualizing results

**Q3**. Create a plot to display the normalized expression of the top 20 genes from the **knockdown versus control** comparison (`res_tableKD`) and **upload the saved imag**e.

**Q4**. For the genes that are differentially expressed in the knockdown versus control comparison (`res_tableKD`), plot an expression heatmap using normalized counts and `pheatmap()` **following the instructions below. Upload the saved image**.

a. The heatmap should only include control and knockdown samples. 

b. Set up a `heat.color`s vector using a palette of your choice from `brewer.pal` (*make sure it is different from the one used in class*).

c. Plot the heatmap without clustering the columns. 

d. Scale expression values by row.

 
## Use significant gene lists to find overlaps between the two comparisons

**Q5**. Using the original cutoff values, perform the following steps:

* padj.cutoff < 0.05

a. Create separate vectors with gene names (Ensembl IDs)  for up-regulated genes and down-regulated genes from `res_tableOE` and save as `up_OE` and `down_O`E, respectively. **Write the code below**:

b. Create separate vectors with gene names (Ensembl IDs) for up-regulated genes and down-regulated genes from `res_tableKD` and save as `up_KD` and `down_KD`, respectively. **Write the code below**:

c. Test for overlaps between the lists:

* How many, and which genes in up_OE are also in down_KD?
* How many, and which genes in up_KD are also in down_OE?
 
