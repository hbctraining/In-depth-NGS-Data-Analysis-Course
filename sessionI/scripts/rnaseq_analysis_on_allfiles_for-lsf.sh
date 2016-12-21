#!/bin/bash

for fq in ~/unix_oct2015/raw_fastq/*.fq
do
  bsub -q priority -n 6 -W 1:30 -R "rusage[mem=4000]" -J rnaseq_mov10 -o %J.out -e %J.err "sh ~/rnaseq_analysis_on_allfiles.sh $fq"
  sleep 1
done
