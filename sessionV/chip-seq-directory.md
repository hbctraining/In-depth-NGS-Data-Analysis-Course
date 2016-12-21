This is the basic directory structure we would like all students to have for chip-seq:

```
~/
├── ngs_course/
    ├── chipseq/
        ├── data/
            ├── trimmed_fastq/
                |-- H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq
            ├── untrimmed_fastq/
            ├── reference_data
                |-- chr12.1.bt2
                |-- chr12.2.bt2
                |-- chr12.3.bt2
                |-- chr12.4.bt2
                |-- chr12.rev.1.bt2
                |-- chr12.rev.2.bt2
            ├── .
        ├── results/
            |-- untrimmed_fastqc
            |-- trimmed_fastqc
            |-- bowtie2
                |-- H1hesc_Input_Rep1_chr12_aln.bam
        ├── .
        └── .
```
