## Session V: ChIP-Seq

### Description

Session V will focus on ChIP-Seq analysis. We will start with read quality control and alignment, then discuss various tools for peak calling, ChIP-Seq QC and visualization for ChIP-Seq analysis. We will talk about the IDR framework when handling replicates with ChIP-seq data and walk-through running IDR on true replicates. Since we have samples from two different IPs we will use DiffBind to demonstrate how to perform differential binding analyses. Finally, we will perform functional analysis of the peak calls, ending our analysis of the ChIP-Seq data. We will end the course by going over a few additional helpful tools. 

**Day1:** We will start the first day with a brief introduction to ChIP-Seq, given by Dr. Shannan Ho Sui, director of the core. Students will then perform quality control and alignment on the sequencing reads. Then, we will perform peak calling using MACS2 followed by QC assessment of the ChIP-Seq data. 

**Day2:** The second day of this session will begin by assesings reproducibility between replicates using bedtools and IDR (Irreproducible Discovery Rate), which will enable us to generate a list of high-confidence peaks. The high-confident peaks will be used for functional analysis using ChipSeeker, the MEME suite and GREAT; these tools help consolidate the output of the analysis with known biology. Next, we will demonstrate the use of DiffBind when comparing peaks between two sample groups in a ChIP-seq experiment and visualize the differentially bound peaks in IGV; this will include a short introduction to bringing in data from Encode for visual comparison. We will end with a discussion of how to integrate ChIP-Seq and RNA-Seq data and explore how to mount O2 on your local machine to efficiently access data on the cluster.   

### Lessons
[Click here for the schedule with links to the lessons.](schedule)

### Learning Objectives
* Understand best practices for designing a ChIP-seq experiment and analysis of the resulting data.
