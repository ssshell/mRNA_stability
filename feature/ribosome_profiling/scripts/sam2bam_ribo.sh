#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/sam
for file in *.sam;do
        base=${file%.*}
        samtools view -h -b -S $file -o ../bam/$base.bam
done
