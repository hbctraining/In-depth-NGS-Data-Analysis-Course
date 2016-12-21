#!/bin/bash

#BSUB -q priority # queue name
#BSUB -W 2:00 # hours:minutes runlimit after which job will be killed.
#BSUB -n 4 # number of cores requested
#BSUB -J rnaseq_mov10_trim         # Job name
#BSUB -o %J.out       # File to which standard out will be written
#BSUB -e %J.err       # File to which standard err will be written

cd ~/unix_oct2015/rnaseq_project/data/untrimmed_fastq

java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
	-threads 4 \
	-phred33 \
	Mov10_oe_1.subset.fq \
	../trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq \
	ILLUMINACLIP:/opt/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 \
	TRAILING:25 \
	MINLEN:35

