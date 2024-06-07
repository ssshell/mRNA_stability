#!/bin/bash

cd ../seq_processing/TotalRibo/
for file in *fastq;do
	base=${file%.*}
	java -jar /Users/serenesun/Documents/HuamingS_WPI/ribosome_profiling_ThisIsTheFinal/seq_processing/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -trimlog ../../scripts/trim_ribo.log $file ../TotalRibo_trim/$base.fastq.gz ILLUMINACLIP:../../scripts/adaptors_SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
done
