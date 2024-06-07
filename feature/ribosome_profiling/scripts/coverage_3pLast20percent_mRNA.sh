#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/bam_filtered
for file in *bam;do
	base=${file%%.*}
	../../../stringtie-2.2.1.OSX_x86_64/stringtie $file -A ../coverage/$base'_'3pLast20percent.tab -G ../../index/msmeg_CombinedAnnotation_CDS_splitByPortions_last20percent.gtf -o ../coverage/$base'_'cds_last20percent_out.gtf
done
