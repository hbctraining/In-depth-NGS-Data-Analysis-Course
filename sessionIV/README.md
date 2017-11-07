### Session IV: RNA-Seq Part III and Related Technologies and Tools

### Description

Session IV starts with a discussion about different methods for differential expression analysis of RNA-seq data. It then segues into a discussion on the considerations for isoform level differential expression (lecture), followed by an interactive lesson on the quantification of isoforms using pseudo-alignment tools and analysis of isoform level differential expression of RNA-seq data. Transcriptome expression abundances for the MOV10 dataset are estimated using [Salmon](https://combine-lab.github.io/salmon/getting_started/). The output of Salmon is used as input to Sleuth to obtain lists of differentially expressed transcripts/genes. The first half of the session ends with a lecture talking about advanced concepts related to using bash, and more specifically using the Orchestra cluster.

The second half of Session IV begins with combining the different components of RNA-seq analysis to create an automated workflow/pipeline using shell scripting. Following this, we introduce HBC's [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/en/latest/) pipeline for RNA-seq analysis with a hands-on demonstration using the Mov10 dataset. Finally, we explore other NGS technologies related to RNA-seq, with a focus on small RNA-seq and single cell RNA-seq (lectures).

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.


### Contents


| Lessons            | Estimated Duration |
|:------------------------|:----------|
| [RNA-seq analysis methods](lectures/) | |
| [Considerations for isoform-level differential expression analysis](lectures/) | |
| [Alignment-free expression estimation using Salmon](lessons/01_salmon.md)| |
| [Isoform-level differential expression using Sleuth](lessons/02_sleuth.md)| |
| [Revisiting Orchestra](lectures/) | |
| [Automating the RNA-seq workflow](lessons/03_automating_workflow.md) | |
| [bcbio-nextgen RNA-seq](lessons/04_bcbio_nextgen.md) | |
| [Accessing genomic data - genome builds and FTP sites](lectures/) | |
| [Other applications of RNA sequencing](lectures/) | |
| [Small RNA-seq analysis](lectures/) | |
| [Single-cell RNA-seq analysis (lecture)](lecture/) and [Seurat RMarkdown demo](lessons/07_single_cell_rnaseq.Rmd) | |



