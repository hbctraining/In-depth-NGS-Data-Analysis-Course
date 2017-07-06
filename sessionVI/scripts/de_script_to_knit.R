## Gene-level differential expression analysis using DESeq2
## June 14th, 2017

## Setup
### Bioconductor and CRAN libraries used
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(reshape)

## Load in data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 
meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
dds <- DESeq(dds)

## QC
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="sampletype")

### Extract the rlog matrix and compute correlation
rld_mat <- assay(rld) 
rld_cor <- cor(rld_mat) 

### Plot heatmap
pheatmap(rld_cor)

### Plot dispersion estimates
plotDispEsts(dds)

## Extract results
res_tableOE <- results(dds, contrast = c("sampletype", "MOV10_overexpression", "control"))
res_tableKD <- results(dds, contrast = c("sampletype", "MOV10_knockdown", "control"))

## Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

## Summarize results
summary(res_tableOE)
summary(res_tableKD)

# Significant results for mov10OE relative to control
sigOE <- subset(data.frame(res_tableOE), padj < padj.cutoff)

# Narrow the results using a log2 foldchange threshold as well
sigOE <- subset(data.frame(res_tableOE), padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Significant results for mov10KD relative to control
sigKD <- subset(data.frame(res_tableKD), padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

## Plot expression for top 20 genes

## Order significant results by padj values
sigOE_ordered <- sigOE[order(sigOE$padj), ]
top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ])

## normalized counts for top 20 significant genes
normalized_counts <- counts(dds, normalized=T)
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

## use melt to modify the format of the data frame
melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))

## add column names that make sense
colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts")

## add metadata to "melted" dataframe
meta$samplename <- rownames(meta)
melted_top20_sigOE <- merge(melted_top20_sigOE, meta)

## plot using ggplot2
ggplot(melted_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
