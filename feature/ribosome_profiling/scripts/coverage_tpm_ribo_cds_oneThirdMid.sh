#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/coverage
for file in *oneThirdMid.tab;do
	base=${file%.*}
	grep -v 'STRG' $file | cut -f1,9| sed 1d |sort > ../coverage_tpm/$base'_'tpm'_'ribo.txt
done
