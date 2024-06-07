#!/bin/bash
cd ../seq_processing/TotalRibo
for file in *.fastq;do
	fastqc $file
	mv *.zip ../fastqc_raw_ribo
	mv *.html ../fastqc_raw_ribo
done
