#!/bin/bash

cd ../seq_processing/TotalmRNA
for file in *.fastq;do
	fastqc $file
	mv *.zip ../fastqc_raw_mrna
	mv *.html ../fastqc_raw_mrna
done
