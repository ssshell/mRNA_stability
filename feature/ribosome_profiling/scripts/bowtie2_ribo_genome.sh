#!/bin/bash

cd ../seq_processing/bowtie2_ribo_tRNArRNA/unaligned/
for file in *fastq;do
	base=${file%.*}
	bowtie2 -x ../../index/msmeg --sensitive-local -p 4 -a -U $file -S ../../bowtie2_ribo_genome/sam/$base.sam
done
