#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/bam_filtered 
for file in *bam;do
	base=${file%%.*}
	../../../stringtie-2.2.1.OSX_x86_64/stringtie $file -A ../coverage/$base'_'oneThird5p.tab -G ../../index/msmeg_CombinedAnnotation_CDS_oneThird_5p.gtf -o ../coverage/$base'_'cds_oneThird5p_out.gtf
done
