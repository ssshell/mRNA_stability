#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/bam_filtered
for file in *bam;do
	base=${file%%.*}
	../../../stringtie-2.2.1.OSX_x86_64/stringtie $file -A ../coverage/$base'_'oneThird3p.tab -G ../../index/msmeg_CombinedAnnotation_CDS_oneThird_3p.gtf -o ../coverage/$base'_'cds_oneThird3p_out.gtf
done
