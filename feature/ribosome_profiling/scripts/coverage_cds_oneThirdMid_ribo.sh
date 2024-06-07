#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/bam_filtered 
for file in *bam;do
	base=${file%%.*}
	../../../stringtie-2.2.1.OSX_x86_64/stringtie $file -A ../coverage/$base'_'oneThirdMid.tab -G ../../index/msmeg_CombinedAnnotation_CDS_oneThird_mid.gtf -o ../coverage/$base'_'cds_oneThirdMid_out.gtf
done
