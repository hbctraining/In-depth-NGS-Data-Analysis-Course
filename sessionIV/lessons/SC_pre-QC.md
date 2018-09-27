# Single-cell RNA-seq: raw sequencing data to counts

Single-cell RNA-seq (scRNA-seq) is an exciting and cutting-edge method for analyzing differences in cellular gene expression, particularly for tissue heterogeneity analyses, lineage tracing, and cell population dynamics. The complexity of scRNA-seq data, which is generally characterized as a large volume of data, representing thousands of cells, and by a low depth of sequencing per cell, resulting in a large number of genes without any corresponding reads, makes analysis of the data more involved than bulk RNA-seq. In addition, the analysis goals can vary depending whether the goal is marker identification, lineage tracing, or some other custom analysis. Therefore, tools specific for scRNA-seq and it's different methods of library preparation are needed. 

The analysis workflow for scRNA-seq is generally similar for the differing scRNA-seq methods, but some specifics regarding the parsing of the UMIs, cell IDs, and sample IDs will differ between them. For example, below is a schematic of the inDrop sequence reads:

<sc_seq_method.png>

While the 10X sequence reads have the UMI and barcodes placed differently:

<10_seq_method.png>

The scRNA-seq method will determine the how to parse the barcodes and UMIs from the sequencing reads. However, the overall workflow will generally follow the same steps regardless of method. The general workflow is shown below:

<sc_workflow.png>

The scRNA-seq method-specific steps are required for the generation of the count matrix, and we will cover what is involved in this later, but after this step, the same methods can be utilized. After generating the count matrix, the raw counts will be assessed to filter out poor quality cells with a low number of genes or UMIs, high mitochondrial gene expression indicative of dying cells, or low number of genes per UMI. After removing the poor quality cells, the cells are clustered based on similarities in transcriptional activity, with the idea that the different cell types separate into the different clusters. After clustering, we can explore genes that are markers for different clusters, which can help identify the cell type of each cluster. Finally, after identification of cell types, there are various types of analyses that can be performed depending on the goal of the experiment.

We are going to start by discussing the first part of this workflow: generating the count matrix from the raw sequencing data.

<sc_gen_matrix_workflow.png>

The sequencing facility will either output the raw sequencing data as BCL format or FASTQ. If the reads are in BCL format, then we will need to convert into FASTQ format. There is a useful tool on O2 called `bcl2fastq` that can easily perform this conversion. We do not demultiplex at this step in the workflow. You may have sequenced 6 samples, but the reads for all samples may be present all in the same BCL or FASTQ file.

<sc_pre-QC_workflow.png>

The generation of the count matrix from the raw sequencing data will go through the following steps for many of the scRNA-seq methods.
The FASTQ files can then be used to parse out the cell barcodes, UMIs, and sample barcodes. Many of the cellular barcodes will match a low number of reads (< 1000 reads) due to encapsulation of free floating RNA from dying cells, small cells, or set of cells that failed for some reason. These excess barcodes need to be filtered out of the sequence data prior to read alignment.


