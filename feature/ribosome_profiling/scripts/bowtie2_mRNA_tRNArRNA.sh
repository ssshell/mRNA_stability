#!/bin/bash

cd ../seq_processing/TotalmRNA
for file in *_1.fastq;do
	base=${file%_*}
	bowtie2 -x ../index/msmeg_tRNArRNA --very-sensitive -p 4 -1 $file -2 $base'_'2.fastq -S ../bowtie2_mRNA_tRNArRNA/sam/$base.sam --un-conc ../bowtie2_mRNA_tRNArRNA/unaligned/$base.fastq
done
