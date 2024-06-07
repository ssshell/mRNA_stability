#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/bam
for file in *.bam;do
	base=${file%.*}
	samtools sort $file -o ../bam_sorted/$base.sorted.bam
done
