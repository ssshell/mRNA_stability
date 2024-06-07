#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/bam_sorted/
for file in *.bam;do
        base=${file%%.*}
        samtools view -b -q 10 -F 1284 $file > ../bam_filtered/$base.filtered.bam
done
