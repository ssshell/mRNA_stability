#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/bam_filtered 
for file in *bam;do
	base=${file%%.*}
	../../../stringtie-2.2.1.OSX_x86_64/stringtie $file -A ../coverage/$base'_'cds_20up.tab -G ../../index/msmeg_CombinedAnnotation_CDS_20up.gtf -o ../coverage/$base'_'cds_20up_out.gtf
done
