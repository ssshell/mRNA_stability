#!/bin/bash

cd ../seq_processing/TotalRibo_trim/
for file in *fastq;do
	base=${file%.*}
	bowtie2 -x ../index/msmeg_tRNArRNA --very-sensitive -p 4 -U $file -S ../bowtie2_ribo_tRNArRNA/sam/$base.sam --un ../bowtie2_ribo_tRNArRNA/unaligned/$base.fastq
done
