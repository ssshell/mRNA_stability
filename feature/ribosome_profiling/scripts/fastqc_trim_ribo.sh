#!/bin/bash
cd ../seq_processing/TotalRibo_trim
for file in *.fastq;do
	fastqc $file
	mv *.zip ../fastqc_trim_ribo
	mv *.html ../fastqc_trim_ribo
done
