#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_tRNArRNA/unaligned/
for file in *.1.fastq;do
	base=${file%%.*}
	bowtie2 -x ../../index/msmeg --sensitive-local -a -p 4 -1 $file -2 $base.2.fastq -S ../../bowtie2_mRNA_genome/sam/$base.sam
done
